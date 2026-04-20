from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Union

import logomaker
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter


PathLike = Union[str, Path]


# ---------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------

def _read_csv(file_path: PathLike) -> pd.DataFrame:
    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"Input file not found: {file_path}")
    return pd.read_csv(file_path)


def _ensure_output_dir(directory: PathLike) -> Path:
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    return directory


def _base_stem(file_path: PathLike) -> str:
    path = Path(file_path)
    stem = path.stem
    if stem.endswith("_processed"):
        stem = stem[:-10]
    return stem


def _save_figure(fig: plt.Figure, directory: PathLike, filename: str) -> Path:
    outdir = _ensure_output_dir(directory)
    outpath = outdir / filename
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)
    return outpath


def _require_columns(df: pd.DataFrame, columns: Sequence[str], context: str = "") -> None:
    missing = [col for col in columns if col not in df.columns]
    if missing:
        prefix = f"{context}: " if context else ""
        raise ValueError(f"{prefix}Missing required columns: {missing}")


def _percentage_columns(df: pd.DataFrame) -> List[str]:
    return [col for col in df.columns if col.startswith("percentage")]


def _filtered_amino_acid_changes(df: pd.DataFrame) -> pd.DataFrame:
    _require_columns(
        df,
        ["Positions within Codon", "Amino Acid", "New_AminoAcid"],
        context="amino-acid change plotting",
    )
    filtered = df.dropna(subset=["Positions within Codon"]).copy()
    filtered["Amino Acid Change"] = (
        filtered["Amino Acid"].astype(str) + "->" + filtered["New_AminoAcid"].astype(str)
    )
    filtered = filtered[filtered["Amino Acid"] != filtered["New_AminoAcid"]]
    return filtered


def _set_theme() -> None:
    sns.set_theme(style="whitegrid")


# ---------------------------------------------------------------------
# Top amino-acid change plot
# ---------------------------------------------------------------------

def create_standard_plot(file_name: PathLike, directory: PathLike, value_threshold: float = 0) -> Path:
    """
    Create a bar chart of the top 10 most frequent amino-acid changes
    above a minimum percentage threshold.
    """
    df = _read_csv(file_name)
    pct_cols = _percentage_columns(df)
    if not pct_cols:
        raise ValueError("No percentage columns found.")

    filtered = df[df[pct_cols].max(axis=1) > value_threshold].copy()
    filtered = filtered.dropna(subset=["Positions within Codon"])
    filtered["Amino Acid Change"] = (
        filtered["Amino Acid"].astype(str) + "→" + filtered["New_AminoAcid"].astype(str)
    )
    filtered = filtered[filtered["Amino Acid"] != filtered["New_AminoAcid"]]

    top_changes = filtered["Amino Acid Change"].value_counts().head(10)
    if top_changes.empty:
        raise ValueError(f"No amino acid changes with percentage > {value_threshold} were found.")

    _set_theme()
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.barplot(
        x=top_changes.index,
        y=top_changes.values,
        palette="viridis",
        edgecolor="black",
        legend=False,
        ax=ax,
    )

    sample_name = _base_stem(file_name)
    ax.set_title(f"Top 10 Most Frequent Amino Acid Changes for {sample_name}", fontsize=16, fontweight="bold")
    ax.set_xlabel("Amino Acid Change", fontsize=14)
    ax.set_ylabel("Frequency", fontsize=14)
    ax.tick_params(axis="x", rotation=45, labelsize=12)
    ax.tick_params(axis="y", labelsize=12)

    for bar, count in zip(ax.patches, top_changes.values):
        ax.annotate(
            f"{count}",
            (bar.get_x() + bar.get_width() / 2, bar.get_height()),
            ha="center",
            va="bottom",
            fontsize=12,
            fontweight="bold",
        )

    fig.tight_layout()
    return _save_figure(fig, directory, f"{sample_name}_barchart_top10AA_Changes.pdf")


# ---------------------------------------------------------------------
# Top mutation positions plot
# ---------------------------------------------------------------------

def create_plot_most_changed_positions(
    file_name: PathLike,
    directory: PathLike,
    value_threshold: float = 0,
) -> List[Path]:
    """
    Create top-20 mutation position plots for each percentage* column.
    """
    df = _read_csv(file_name)
    _require_columns(
        df,
        ["effect", "position", "Amino Acid", "New_AminoAcid", "Gene"],
        context="create_plot_most_changed_positions",
    )

    pct_cols = _percentage_columns(df)
    if not pct_cols:
        raise ValueError("No percentage columns found.")

    effect_df = df[
        (df[pct_cols] > value_threshold).any(axis=1)
        & df["effect"].notna()
        & (df["effect"] != "")
    ].copy()

    _set_theme()
    color_palette = {"AA_Change": "#D55E00", "silent": "#009E73"}
    outputs: List[Path] = []
    sample_name = _base_stem(file_name)

    for col in pct_cols:
        top_20 = effect_df.nlargest(20, col).sort_values("position")
        if top_20.empty:
            continue

        positions = top_20["position"].astype(str)
        values = top_20[col]
        colors = top_20["effect"].map(color_palette).fillna("#999999")

        fig, ax = plt.subplots(figsize=(12, 8))
        bars = ax.bar(range(len(positions)), values, color=colors, edgecolor="black", linewidth=0.7)

        ax.set_title(f"Top 20 Mutation Positions by Percentage for {sample_name}", fontsize=18, fontweight="bold")
        ax.set_xlabel("Position", fontsize=14, labelpad=10)
        ax.set_ylabel("Mutation Percentage (%)", fontsize=14)
        ax.set_xticks(range(len(positions)))
        ax.set_xticklabels(positions, rotation=45, ha="center", fontsize=11)
        ax.tick_params(axis="y", labelsize=12)
        ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: int(x)))

        for bar, amino_acid, new_amino_acid, effect in zip(
            bars,
            top_20["Amino Acid"],
            top_20["New_AminoAcid"],
            top_20["effect"],
        ):
            if effect == "AA_Change":
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.5,
                    f"{amino_acid}>{new_amino_acid}",
                    ha="center",
                    va="bottom",
                    fontsize=10,
                    color="black",
                )

        if "PositionProtein" in top_20.columns:
            protein_labels = top_20.apply(
                lambda row: (
                    f"{row['Gene']} ({row['Amino Acid']}"
                    f"{int(row['PositionProtein']) if pd.notna(row['PositionProtein']) else 'N/A'}"
                    f"{row['New_AminoAcid']})"
                    if row["effect"] == "AA_Change"
                    else row["Gene"]
                ),
                axis=1,
            )
            ax2 = ax.secondary_xaxis("top")
            ax2.set_ticks(range(len(positions)), labels=protein_labels, rotation=45, ha="left", fontsize=8)
            ax2.set_xlim(ax.get_xlim())
            ax2.spines["top"].set_linewidth(1.2)
            ax2.set_xlabel("")

        aa_change_patch = mpatches.Patch(color=color_palette["AA_Change"], label="AA Change")
        silent_patch = mpatches.Patch(color=color_palette["silent"], label="Silent Mutation")
        ax.legend(handles=[aa_change_patch, silent_patch], title="Mutation Type", loc="upper right", frameon=True)

        fig.tight_layout()
        outputs.append(_save_figure(fig, directory, f"{sample_name}_{col}_barchart_top20_Mutations.pdf"))

    return outputs


# ---------------------------------------------------------------------
# Pie and donut charts
# ---------------------------------------------------------------------

def create_pie_chart(filename: PathLike, directory: PathLike, value_threshold: float = 0) -> Path:
    _read = _read_csv(filename)
    _require_columns(_read, ["Positions within Codon", "percentage"], context="create_pie_chart")

    valid_rows = _read[
        (_read["percentage"] >= value_threshold)
        & _read["Positions within Codon"].notna()
    ]
    if valid_rows.empty:
        raise ValueError(f"No valid rows after filtering in {filename}")

    position_counts = valid_rows["Positions within Codon"].value_counts().sort_index()
    labels = [f"Position {int(pos)}" for pos in position_counts.index]

    def autopct_func(pct: float, allvals: pd.Series) -> str:
        absolute = int(np.round(pct / 100.0 * np.sum(allvals)))
        return f"{pct:.1f}% ({absolute})"

    _set_theme()
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.pie(
        position_counts,
        labels=labels,
        autopct=lambda pct: autopct_func(pct, position_counts),
        startangle=140,
        colors=sns.color_palette("pastel")[0:3],
        wedgeprops={"edgecolor": "black", "linewidth": 0.8},
    )
    ax.set_title(
        f"Distribution of Codon Positions\nfor {_base_stem(filename)} (≥ {value_threshold}% threshold)",
        fontsize=14,
        fontweight="bold",
    )

    return _save_figure(fig, directory, f"{_base_stem(filename)}_positions_pie_chart.pdf")


def create_donut_chart(filename: PathLike, directory: PathLike, value_threshold: float = 0) -> Path:
    df = _read_csv(filename)
    _require_columns(df, ["Positions within Codon", "percentage"], context="create_donut_chart")

    valid_rows = df[
        (df["percentage"] >= value_threshold)
        & df["Positions within Codon"].notna()
    ]
    if valid_rows.empty:
        raise ValueError(f"No valid rows after filtering in {filename}")

    position_counts = valid_rows["Positions within Codon"].value_counts().sort_index()
    labels = [f"Position {int(pos)}" for pos in position_counts.index]

    def autopct_func(pct: float, allvals: pd.Series) -> str:
        absolute = int(np.round(pct / 100.0 * np.sum(allvals)))
        return f"{pct:.1f}%\n({absolute})"

    _set_theme()
    fig, ax = plt.subplots(figsize=(8, 8))
    wedges, texts, autotexts = ax.pie(
        position_counts,
        labels=labels,
        autopct=lambda pct: autopct_func(pct, position_counts),
        startangle=140,
        colors=sns.color_palette("pastel")[0:3],
        wedgeprops={"edgecolor": "black", "linewidth": 0.8},
        textprops={"fontsize": 12, "weight": "bold"},
    )

    center_circle = plt.Circle((0, 0), 0.70, facecolor="white", edgecolor="black", linewidth=0.8)
    ax.add_artist(center_circle)

    sample_name = _base_stem(filename)
    ax.text(0, 0, sample_name, ha="center", va="center", fontsize=16, fontweight="bold")

    for text in autotexts:
        text.set_fontsize(11)
        text.set_color("black")
        text.set_fontweight("normal")
        text.set_verticalalignment("top")

    ax.set_title(f"Distribution of Codon Positions\n(≥ {value_threshold}% threshold)", fontsize=14, fontweight="bold")
    fig.tight_layout()

    return _save_figure(fig, directory, f"{_base_stem(filename)}_positions_donut_chart.pdf")


# ---------------------------------------------------------------------
# Range-based amino-acid change plot
# ---------------------------------------------------------------------

def create_ranges_plot(file_name: PathLike, directory: PathLike, percentage_ranges: List[Dict[str, object]]) -> Path:
    df = _read_csv(file_name)
    filtered = _filtered_amino_acid_changes(df)

    _set_theme()
    max_bars = 20
    bar_height = 0.6

    subplot_heights = [
        max(
            len(filtered[rng["condition"](filtered)]["Amino Acid Change"].value_counts().head(max_bars)),
            1,
        ) * bar_height
        for rng in percentage_ranges
    ]
    total_height = sum(subplot_heights) + len(percentage_ranges)

    fig, axes = plt.subplots(
        len(percentage_ranges),
        1,
        figsize=(12, total_height),
        gridspec_kw={"height_ratios": subplot_heights},
    )

    if len(percentage_ranges) == 1:
        axes = [axes]

    for idx, rng in enumerate(percentage_ranges):
        range_df = filtered[rng["condition"](filtered)]
        aa_change_counts = range_df["Amino Acid Change"].value_counts().head(max_bars)

        sns.barplot(
            x=aa_change_counts.values,
            y=aa_change_counts.index,
            ax=axes[idx],
            palette=["#D55E00"] * len(aa_change_counts),
            dodge=False,
        )

        axes[idx].set_ylim(-0.5, max(len(aa_change_counts), 1) - 0.5)
        axes[idx].set_title(rng["title"], fontsize=14, fontweight="bold")
        axes[idx].set_xlabel("Frequency", fontsize=12)
        axes[idx].set_ylabel("Amino Acid Change", fontsize=12)

    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.suptitle(f"Top 20 Amino Acid Changes for {_base_stem(file_name)}", fontsize=16, fontweight="bold")

    return _save_figure(fig, directory, f"{_base_stem(file_name)}_top_positions_ranges_chart.pdf")


# ---------------------------------------------------------------------
# Bubble plots
# ---------------------------------------------------------------------

def create_bubble_plot(filename: PathLike, directory: PathLike) -> Path:
    df = _read_csv(filename)
    filtered = _filtered_amino_acid_changes(df)

    change_counts = filtered.groupby(["Amino Acid", "New_AminoAcid"]).size().reset_index(name="count")
    amino_acids = sorted(filtered["Amino Acid"].unique())

    fig, ax = plt.subplots(figsize=(12, 8))
    ax.scatter(
        change_counts["New_AminoAcid"],
        change_counts["Amino Acid"],
        s=change_counts["count"] * 10,
        alpha=0.5,
    )

    ax.set_xlabel("New Amino Acid")
    ax.set_ylabel("Amino Acid")
    ax.set_title(f"Amino Acid Changes for {_base_stem(filename)}")
    ax.grid(True)

    ax.set_xticks(range(len(amino_acids)))
    ax.set_xticklabels(amino_acids)
    ax.set_yticks(range(len(amino_acids)))
    ax.set_yticklabels(amino_acids)

    max_count = int(change_counts["count"].max())
    sizes = [1, max(1, max_count // 2), max_count]
    legend_bubbles = [
        ax.scatter([], [], s=size * 10, alpha=0.5, color="white", edgecolor="black", label=f"{size} occurrences")
        for size in sizes
    ]
    ax.legend(handles=legend_bubbles, scatterpoints=1, frameon=True, labelspacing=1, title="Bubble Size")

    return _save_figure(fig, directory, f"{_base_stem(filename)}_amino_acid_changes_bubble_plot.pdf")


def create_bubble_ranges_plot(
    filename: PathLike,
    directory: PathLike,
    percentage_ranges: List[Dict[str, object]],
) -> List[Path]:
    df = _read_csv(filename)
    filtered = _filtered_amino_acid_changes(df)
    amino_acids = sorted(filtered["Amino Acid"].unique())

    outputs: List[Path] = []

    for rng in percentage_ranges:
        range_df = filtered[rng["condition"](filtered)]
        change_counts = range_df.groupby(["Amino Acid", "New_AminoAcid"]).size().reset_index(name="count")
        if change_counts.empty:
            continue

        fig, ax = plt.subplots(figsize=(12, 8))
        ax.scatter(
            change_counts["New_AminoAcid"],
            change_counts["Amino Acid"],
            s=change_counts["count"] * 10,
            alpha=0.5,
        )

        base_filename = _base_stem(filename)
        ax.set_xlabel("New Amino Acid")
        ax.set_ylabel("Amino Acid")
        ax.set_title(f'{rng["title"]} for {base_filename}')
        ax.grid(True)

        ax.set_xticks(range(len(amino_acids)))
        ax.set_xticklabels(amino_acids)
        ax.set_yticks(range(len(amino_acids)))
        ax.set_yticklabels(amino_acids)

        max_count = int(change_counts["count"].max())
        sizes = [1, max(1, max_count // 2), max_count]
        legend_bubbles = [
            ax.scatter([], [], s=size * 10, alpha=0.5, color="white", edgecolor="black", label=f"{size} occurrences")
            for size in sizes
        ]
        ax.legend(handles=legend_bubbles, scatterpoints=1, frameon=True, labelspacing=1, title="Bubble Size")

        outputs.append(
            _save_figure(
                fig,
                directory,
                f'{base_filename}_{rng["filename_suffix"]}_amino_acid_changes_bubble_plot.pdf',
            )
        )

    return outputs


def create_bubble_ranges_1fig_plot(
    filename: PathLike,
    directory: PathLike,
    percentage_ranges: List[Dict[str, object]],
) -> Path:
    df = _read_csv(filename)
    filtered = _filtered_amino_acid_changes(df)
    amino_acids = sorted(filtered["Amino Acid"].unique())

    fig, axs = plt.subplots(1, len(percentage_ranges), figsize=(18, 6), sharey=True)
    if len(percentage_ranges) == 1:
        axs = [axs]

    for i, rng in enumerate(percentage_ranges):
        range_df = filtered[rng["condition"](filtered)]
        change_counts = range_df.groupby(["Amino Acid", "New_AminoAcid"]).size().reset_index(name="count")

        axs[i].scatter(
            change_counts["New_AminoAcid"],
            change_counts["Amino Acid"],
            s=change_counts["count"] * 10,
            alpha=0.5,
        )

        axs[i].set_xlabel("New Amino Acid")
        axs[i].set_title(rng["title"])
        axs[i].grid(True)
        axs[i].set_xticks(range(len(amino_acids)))
        axs[i].set_xticklabels(amino_acids, rotation=90)
        axs[i].set_yticks(range(len(amino_acids)))
        axs[i].set_yticklabels(amino_acids)
        if i == 0:
            axs[i].set_ylabel("Amino Acid")

    max_count = int(filtered.groupby(["Amino Acid", "New_AminoAcid"]).size().max())
    sizes = [1, max(1, max_count // 2), max_count]
    legend_bubbles = [
        plt.scatter([], [], s=size * 10, alpha=0.5, color="white", edgecolor="black", label=f"{size} occurrences")
        for size in sizes
    ]
    fig.legend(handles=legend_bubbles, scatterpoints=1, frameon=True, labelspacing=1, title="Bubble Size", loc="upper right")

    plt.tight_layout(rect=[0, 0, 0.9, 0.95])
    plt.suptitle(f"Amino Acid Changes Frequency for {_base_stem(filename)}", y=1.02)

    return _save_figure(fig, directory, f"{_base_stem(filename)}_combined_amino_acid_changes_bubble_plot_2.pdf")


def create_combined_bubble_plot(
    filename: PathLike,
    directory: PathLike,
    percentage_ranges: List[Dict[str, object]],
) -> Path:
    df = _read_csv(filename)
    filtered = _filtered_amino_acid_changes(df)
    amino_acids = sorted(filtered["Amino Acid"].unique())

    fig, ax = plt.subplots(figsize=(12, 8))

    for rng in percentage_ranges:
        range_df = filtered[rng["condition"]]
        change_counts = range_df.groupby(["Amino Acid", "New_AminoAcid"]).size().reset_index(name="count")

        ax.scatter(
            change_counts["New_AminoAcid"],
            change_counts["Amino Acid"],
            s=change_counts["count"] * 10,
            alpha=0.5,
            label=rng["title"],
            color=rng["color"],
        )

    ax.set_xlabel("New Amino Acid")
    ax.set_ylabel("Amino Acid")
    ax.set_title(f"Amino Acid Changes Frequency for {_base_stem(filename)}")
    ax.grid(True)

    ax.set_xticks(range(len(amino_acids)))
    ax.set_xticklabels(amino_acids, rotation=90)
    ax.set_yticks(range(len(amino_acids)))
    ax.set_yticklabels(amino_acids)

    range_legend = ax.legend(title="Percentage Range", bbox_to_anchor=(1.05, 1), loc="upper left", frameon=True)

    max_count = int(filtered.groupby(["Amino Acid", "New_AminoAcid"]).size().max())
    sizes = [1, max(1, max_count // 2), max_count]
    legend_bubbles = [
        plt.scatter([], [], s=size * 10, alpha=0.5, color="white", edgecolor="black", label=f"{size} occurrences")
        for size in sizes
    ]

    fig.legend(
        handles=legend_bubbles,
        scatterpoints=1,
        frameon=True,
        labelspacing=1,
        title="Bubble Size",
        bbox_to_anchor=(1.05, 0.75),
        loc="upper left",
    )
    ax.add_artist(range_legend)

    plt.tight_layout(rect=[0, 0, 0.85, 1])
    return _save_figure(fig, directory, f"{_base_stem(filename)}_combined_amino_acid_changes_bubble_plot.pdf")


# ---------------------------------------------------------------------
# Overview plots
# ---------------------------------------------------------------------

def create_overview_plot_counts(
    results_df: pd.DataFrame,
    directory: PathLike,
    percentage_ranges: List[Dict[str, object]],
) -> Path:
    effect_types = ["Silent Count", "AA Change Count"]
    colors = {"AA_Change": "#D55E00", "silent": "#009E73"}

    max_y = results_df[effect_types].max().max()
    fig, axes = plt.subplots(1, len(percentage_ranges), figsize=(24, 6))
    if len(percentage_ranges) == 1:
        axes = [axes]

    for i, rng in enumerate(percentage_ranges):
        range_title = rng["title"]
        subset = results_df[results_df["Range"] == range_title].sort_values(by="File")

        if not subset.empty:
            x = np.arange(len(subset["File"]))
            width = 0.35

            ax = axes[i]
            ax.bar(x - width / 2, subset["Silent Count"], width, label="Silent", color=colors["silent"])
            ax.bar(x + width / 2, subset["AA Change Count"], width, label="AA Change", color=colors["AA_Change"])

            ax.set_title(range_title)
            ax.set_ylabel("Count")
            ax.set_xticks(x)
            ax.set_xticklabels(subset["File"], rotation=90)
            ax.set_ylim(0, max_y)

            if i == len(percentage_ranges) - 1:
                ax.legend(loc="upper left")

    plt.tight_layout()
    return _save_figure(fig, directory, "effect_counts_plot.pdf")


def create_overview_plot_percentages(
    results_df: pd.DataFrame,
    directory: PathLike,
    percentage_ranges: List[Dict[str, object]],
) -> Path:
    effect_types = ["Silent Percentage", "AA Change Percentage"]
    colors = {"AA_Change": "#D55E00", "silent": "#009E73"}

    max_y = results_df[effect_types].max().max()
    fig, axes = plt.subplots(1, len(percentage_ranges), figsize=(24, 6))
    if len(percentage_ranges) == 1:
        axes = [axes]

    for i, rng in enumerate(percentage_ranges):
        range_title = rng["title"]
        subset = results_df[results_df["Range"] == range_title].sort_values(by="File")

        if not subset.empty:
            x = np.arange(len(subset["File"]))
            width = 0.35

            ax = axes[i]
            ax.bar(x - width / 2, subset["Silent Percentage"], width, label="Silent Percentage", color=colors["silent"])
            ax.bar(x + width / 2, subset["AA Change Percentage"], width, label="AA Change Percentage", color=colors["AA_Change"])

            ax.set_title(range_title)
            ax.set_ylabel("Percentage")
            ax.set_xticks(x)
            ax.set_xticklabels(subset["File"], rotation=90)
            ax.set_ylim(0, max_y)

            if i == len(percentage_ranges) - 1:
                ax.legend(loc="upper left")

    plt.tight_layout()
    return _save_figure(fig, directory, "effect_percentages_plot.pdf")


# ---------------------------------------------------------------------
# Percentage distribution plots
# ---------------------------------------------------------------------

def create_histogram_per_sample(
    filename: PathLike,
    directory: PathLike,
    min_percentage: int = 0,
    max_percentage: int = 100,
) -> Path:
    df = _read_csv(filename)
    if "percentage" not in df.columns:
        raise ValueError(f"Skipping {filename}, no 'percentage' column found.")

    filtered = df.dropna(subset=["percentage"])
    filtered = filtered[
        (filtered["percentage"] >= min_percentage)
        & (filtered["percentage"] <= max_percentage)
    ]

    if filtered.empty:
        raise ValueError(f"No valid percentage data in {filename}.")

    fig, ax = plt.subplots(figsize=(10, 6))
    bins = list(range(min_percentage, max_percentage + 1, 1))
    ax.hist(filtered["percentage"], bins=bins, edgecolor="black", alpha=0.7)

    ax.set_title(f"Percentage Distribution ({min_percentage}-{max_percentage}%) for {_base_stem(filename)}")
    ax.set_xlabel("Percentage")
    ax.set_ylabel("Count")

    return _save_figure(fig, directory, f"{_base_stem(filename)}_percentage_histogram.pdf")


def create_line_plot_per_sample(
    filename: PathLike,
    directory: PathLike,
    min_percentage: int = 0,
    max_percentage: int = 100,
    bin_width: int = 1,
) -> Path:
    df = _read_csv(filename)
    if "percentage" not in df.columns:
        raise ValueError(f"Skipping {filename}, no 'percentage' column found.")

    filtered = df.dropna(subset=["percentage"])
    filtered = filtered[
        (filtered["percentage"] >= min_percentage)
        & (filtered["percentage"] <= max_percentage)
    ]

    if filtered.empty:
        raise ValueError(
            f"No valid percentage data in the range {min_percentage}-{max_percentage}% for {filename}."
        )

    bins = np.arange(min_percentage, max_percentage + 1, bin_width)
    counts, edges = np.histogram(filtered["percentage"], bins=bins)
    bin_centers = (edges[:-1] + edges[1:]) / 2

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(bin_centers, counts, marker="o", markersize=2, linestyle="-", color="b", alpha=0.7)

    ax.set_title(f"Percentage Distribution ({min_percentage}-{max_percentage}%) for {_base_stem(filename)}")
    ax.set_xlabel("Percentage")
    ax.set_ylabel("Count")

    return _save_figure(
        fig,
        directory,
        f"{_base_stem(filename)}_{min_percentage}_to_{max_percentage}_percentage_line_plot.pdf",
    )


# ---------------------------------------------------------------------
# All RNA types plot
# ---------------------------------------------------------------------

def create_plot_all_rna_types(file_name: PathLike, directory: PathLike) -> List[Path]:
    df = _read_csv(file_name)
    _require_columns(
        df,
        ["position", "RNA Type", "effect", "Gene", "Amino Acid", "New_AminoAcid"],
        context="create_plot_all_rna_types",
    )

    pct_cols = _percentage_columns(df)
    if not pct_cols:
        raise ValueError("No percentage columns found.")

    _set_theme()
    color_palette = {
        "AA_Change_mRNA": "#D55E00",
        "silent_mRNA": "#009E73",
        "tRNA": "#0072B2",
        "rRNA": "#CC79A7",
        "other": "#000000",
    }

    outputs: List[Path] = []
    sample_name = _base_stem(file_name)

    for col in pct_cols:
        top_20 = df.nlargest(20, col).sort_values("position")
        if top_20.empty:
            continue

        positions = top_20["position"].astype(str)
        values = top_20[col]

        fig, ax = plt.subplots(figsize=(12, 8))
        colors = top_20.apply(
            lambda row: (
                color_palette["AA_Change_mRNA"] if row["RNA Type"] == "mRNA" and row["effect"] == "AA_Change"
                else color_palette["silent_mRNA"] if row["RNA Type"] == "mRNA" and row["effect"] == "silent"
                else color_palette["tRNA"] if row["RNA Type"] == "tRNA"
                else color_palette["rRNA"] if row["RNA Type"] == "rRNA"
                else color_palette["other"]
            ),
            axis=1,
        )

        bars = ax.bar(range(len(positions)), values, color=colors, edgecolor="black", linewidth=0.7)

        ax.set_title(f"Top 20 Mutation positions for {sample_name}", fontsize=18, fontweight="bold")
        ax.set_xlabel("Position", fontsize=14, labelpad=10)
        ax.set_ylabel("Mutation Percentage (%)", fontsize=14)
        ax.set_xticks(range(len(positions)))
        ax.set_xticklabels(positions, rotation=45, ha="center", fontsize=11)
        ax.tick_params(axis="y", labelsize=12)
        ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: int(x)))

        for bar, rna_type, amino_acid, new_amino_acid in zip(
            bars,
            top_20["RNA Type"],
            top_20["Amino Acid"],
            top_20["New_AminoAcid"],
        ):
            if rna_type == "mRNA":
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.5,
                    f"{amino_acid}>{new_amino_acid}",
                    ha="center",
                    va="bottom",
                    fontsize=10,
                    color="black",
                )

        if "PositionProtein" in top_20.columns:
            protein_labels = top_20.apply(
                lambda row: (
                    f"{row['Gene']} ({row['Amino Acid']}"
                    f"{int(row['PositionProtein']) if pd.notna(row['PositionProtein']) else 'N/A'}"
                    f"{row['New_AminoAcid']})"
                    if row["RNA Type"] == "mRNA"
                    else f"{row['Gene']}"
                ),
                axis=1,
            )
            ax2 = ax.secondary_xaxis("top")
            ax2.set_ticks(range(len(positions)), labels=protein_labels, rotation=45, ha="left", fontsize=8)
            ax2.set_xlim(ax.get_xlim())
            ax2.spines["top"].set_linewidth(1.2)
            ax2.set_xlabel("")

        aa_change_patch = mpatches.Patch(color=color_palette["AA_Change_mRNA"], label="AA Change (mRNA)")
        silent_patch = mpatches.Patch(color=color_palette["silent_mRNA"], label="Silent (mRNA)")
        trna_patch = mpatches.Patch(color=color_palette["tRNA"], label="tRNA")
        rrna_patch = mpatches.Patch(color=color_palette["rRNA"], label="rRNA")
        ax.legend(
            handles=[aa_change_patch, silent_patch, trna_patch, rrna_patch],
            title="Mutation Type",
            loc="upper right",
            frameon=True,
        )

        fig.tight_layout()
        outputs.append(_save_figure(fig, directory, f"{sample_name}_{col}_barchart_all_positions.pdf"))

    return outputs


# ---------------------------------------------------------------------
# Sequence logo
# ---------------------------------------------------------------------

def create_sequence_logo_from_percentage(
    file_name: PathLike,
    directory: PathLike,
    value_threshold: float,
) -> Path:
    df = _read_csv(file_name)
    _require_columns(
        df,
        ["percentage", "position", "Original Nucleotide"],
        context="create_sequence_logo_from_percentage",
    )

    high_percentage_rows = df[df["percentage"] > value_threshold]

    filtered_sequences: List[str] = []
    for _, row in high_percentage_rows.iterrows():
        pos = int(row["position"]) - 1
        surrounding_sequence = df.iloc[max(0, pos - 5): pos + 6]["Original Nucleotide"].astype(str).str.cat()

        if len(surrounding_sequence) == 11:
            filtered_sequences.append(surrounding_sequence)

    if not filtered_sequences:
        raise ValueError("No sequences found with the given criteria.")

    count_matrix = {nucleotide: [0] * 11 for nucleotide in "ATGC"}
    for sequence in filtered_sequences:
        for i, nucleotide in enumerate(sequence):
            if nucleotide in count_matrix:
                count_matrix[nucleotide][i] += 1

    count_df = pd.DataFrame(count_matrix)
    frequency_df = count_df.div(count_df.sum(axis=1), axis=0)

    fig, ax = plt.subplots(figsize=(10, 6))
    logo = logomaker.Logo(frequency_df, ax=ax, shade_below=0.5, fade_below=0.5)
    logo.style_xticks(rotation=45, fmt="%d", anchor=0)
    logo.ax.set_title(f"Sequence Logo for High Percentage Sequences (> {value_threshold}%)")
    logo.ax.set_ylabel("Frequency")

    return _save_figure(fig, directory, f"{_base_stem(file_name)}_sequence_logo.pdf")


# ---------------------------------------------------------------------
# Gene-specific plots
# ---------------------------------------------------------------------

def create_fig_gene_mut(
    file_path: PathLike,
    directory: PathLike,
    gene: str,
    start_offset: int = 10,
    end_offset: int = 10,
    highlight: Optional[Union[int, Tuple[int, int]]] = None,
) -> Path:
    df = _read_csv(file_path)
    _require_columns(df, ["Gene", "position", "percentage"], context="create_fig_gene_mut")

    df = df.copy()
    df["Gene"] = df["Gene"].astype(str).str.strip().str.upper()
    gene = gene.strip().upper()

    gene_data = df[df["Gene"] == gene]
    if gene_data.empty:
        raise ValueError(f"No data found for gene: {gene}")

    min_position = gene_data["position"].min()
    max_position = gene_data["position"].max()

    start_position = max(0, min_position - start_offset)
    end_position = max_position + end_offset

    offset_data = df[(df["position"] >= start_position) & (df["position"] <= end_position)]
    x_values = offset_data["position"]
    y_values = offset_data["percentage"]

    range_length = end_position - start_position
    interval = 10 if range_length < 200 else 100

    fig, ax = plt.subplots(figsize=(14, 7))

    ax.axvspan(min_position, max_position, color="lightblue", alpha=0.3, label="Gene Region")

    if highlight:
        highlight_start, highlight_end = highlight if isinstance(highlight, tuple) else (highlight, highlight)
        ax.axvspan(highlight_start, highlight_end, color="yellow", alpha=0.3, label="Highlighted Region")

    ax.bar(x_values, y_values, color="lightsteelblue", edgecolor="black", linewidth=0.7)

    x_tick_start = (start_position // interval) * interval
    xticks = [pos for pos in range(x_tick_start, end_position + 1, interval)]
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticks, rotation=45, ha="right", fontsize=10)
    ax.set_xlim(start_position, end_position)

    ax.set_xlabel("Position", fontsize=12, labelpad=10)
    ax.set_ylabel("Percentage (%)", fontsize=12, labelpad=10)
    ax.set_ylim(0, y_values.max() * 1.2)
    ax.grid(axis="y", linestyle="--", linewidth=0.5, alpha=0.7)

    base_name = _base_stem(file_path)
    ax.set_title(base_name, fontsize=16, fontweight="bold", pad=20)

    ax.text(
        x=(min_position + max_position) / 2,
        y=y_values.max() * 1.1,
        s=f"Gene: {gene}",
        ha="center",
        va="center",
        fontsize=14,
        fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.3", edgecolor="black", facecolor="lightgrey"),
    )

    fig.tight_layout()
    return _save_figure(fig, directory, f"{base_name}_{gene}_mutation_percentage_barchart.pdf")


def create_fig_gene_codon(
    file_path: PathLike,
    directory: PathLike,
    gene: str,
    start_offset: int = 10,
    end_offset: int = 10,
    min_threshold: float = 0,
    highlight: Optional[Union[int, Tuple[int, int]]] = None,
) -> Path:
    df = _read_csv(file_path)
    _require_columns(
        df,
        ["Gene", "position", "percentage", "effect", "RNA Type", "Amino Acid", "New_AminoAcid"],
        context="create_fig_gene_codon",
    )

    df = df.copy()
    df["Gene"] = df["Gene"].astype(str).str.strip().str.upper()
    gene = gene.strip().upper()

    gene_data = df[df["Gene"] == gene]
    if gene_data.empty:
        raise ValueError(f"No data found for gene: {gene}")

    gene_data = gene_data[gene_data["percentage"] >= min_threshold]
    if gene_data.empty:
        raise ValueError(f"No data above the minimum threshold for gene: {gene}")

    min_position = gene_data["position"].min()
    max_position = gene_data["position"].max()

    start_position = max(0, min_position - start_offset)
    end_position = max_position + end_offset

    offset_data = df[(df["position"] >= start_position) & (df["position"] <= end_position)]
    offset_data = offset_data[offset_data["percentage"] >= min_threshold]

    x_values = offset_data["position"]
    y_values = offset_data["percentage"]

    def get_bar_color(row: pd.Series) -> str:
        if row["RNA Type"] == "mRNA":
            return "orange" if row["effect"] == "AA_Change" else "green"
        if row["RNA Type"] == "tRNA":
            return "blue"
        if row["RNA Type"] == "rRNA":
            return "pink"
        return "grey"

    colors = offset_data.apply(get_bar_color, axis=1)

    range_length = end_position - start_position
    interval = 10 if range_length < 200 else 100

    fig, ax = plt.subplots(figsize=(14, 7))
    ax.axvspan(min_position, max_position, color="lightblue", alpha=0.3, label="Gene Region")

    if highlight:
        highlight_start, highlight_end = highlight if isinstance(highlight, tuple) else (highlight, highlight)
        ax.axvspan(highlight_start, highlight_end, color="yellow", alpha=0.3, label="Highlighted Region")

    bars = ax.bar(x_values, y_values, color=colors, edgecolor="black", linewidth=0.7)

    if "PositionProtein" in offset_data.columns:
        for bar, effect, amino_acid, new_amino_acid, position_protein in zip(
            bars,
            offset_data["effect"],
            offset_data["Amino Acid"],
            offset_data["New_AminoAcid"],
            offset_data["PositionProtein"],
        ):
            if effect == "AA_Change" and pd.notna(position_protein):
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.1,
                    f"p.{amino_acid}{int(position_protein)}{new_amino_acid}",
                    ha="center",
                    va="bottom",
                    fontsize=10,
                    color="black",
                )

    x_tick_start = (start_position // interval) * interval
    xticks = [pos for pos in range(x_tick_start, end_position + 1, interval)]
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticks, rotation=45, ha="right", fontsize=10)
    ax.set_xlim(start_position, end_position)

    ax.set_xlabel("Position", fontsize=12, labelpad=10)
    ax.set_ylabel("Percentage (%)", fontsize=12, labelpad=10)
    ax.set_ylim(0, y_values.max() * 1.2)
    ax.grid(axis="y", linestyle="--", linewidth=0.5, alpha=0.7)

    base_name = _base_stem(file_path)
    ax.set_title(base_name, fontsize=16, fontweight="bold", pad=20)

    ax.text(
        x=(min_position + max_position) / 2,
        y=y_values.max() * 1.1,
        s=f"Gene: {gene}",
        ha="center",
        va="center",
        fontsize=14,
        fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.3", edgecolor="black", facecolor="lightgrey"),
    )

    legend_handles = [
        Patch(facecolor="orange", edgecolor="black", linewidth=0.7, label="mRNA: AA Change"),
        Patch(facecolor="green", edgecolor="black", linewidth=0.7, label="mRNA: Silent"),
        Patch(facecolor="blue", edgecolor="black", linewidth=0.7, label="tRNA"),
        Patch(facecolor="pink", edgecolor="black", linewidth=0.7, label="rRNA"),
    ]
    ax.legend(handles=legend_handles, title="Bar Colors", loc="upper right", frameon=True)

    fig.tight_layout()
    return _save_figure(fig, directory, f"{base_name}_{gene}_codon_percentage_barchart.pdf")


# ---------------------------------------------------------------------
# Backward-compatible aliases
# ---------------------------------------------------------------------

create_Standard_plot = create_standard_plot
create_plot_MostChangedPositions = create_plot_most_changed_positions
create_plot_AllRNA_Types = create_plot_all_rna_types


__all__ = [
    "create_standard_plot",
    "create_plot_most_changed_positions",
    "create_pie_chart",
    "create_donut_chart",
    "create_ranges_plot",
    "create_bubble_plot",
    "create_bubble_ranges_plot",
    "create_bubble_ranges_1fig_plot",
    "create_combined_bubble_plot",
    "create_overview_plot_counts",
    "create_overview_plot_percentages",
    "create_histogram_per_sample",
    "create_line_plot_per_sample",
    "create_plot_all_rna_types",
    "create_sequence_logo_from_percentage",
    "create_fig_gene_mut",
    "create_fig_gene_codon",
    "create_Standard_plot",
    "create_plot_MostChangedPositions",
    "create_plot_AllRNA_Types",
]
