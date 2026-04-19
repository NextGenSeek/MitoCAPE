from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Sequence, Union

import pandas as pd

from .references import build_reference_dataframe, get_reference_data
from .processing import (
    add_gnomad_metrics,
    calculate_effect_counts_and_percentages,
    calculate_effect_counts_and_percentages_ranges,
    process_file,
)
from .plotting import (
    create_standard_plot,
    create_plot_most_changed_positions,
    create_pie_chart,
    create_donut_chart,
    create_ranges_plot,
    create_bubble_plot,
    create_bubble_ranges_plot,
    create_bubble_ranges_1fig_plot,
    create_combined_bubble_plot,
    create_overview_plot_counts,
    create_overview_plot_percentages,
    create_histogram_per_sample,
    create_line_plot_per_sample,
    create_plot_all_rna_types,
    create_sequence_logo_from_percentage,
    create_fig_gene_mut,
    create_fig_gene_codon,
)


PathLike = Union[str, Path]


def run_analysis(
    input_file: PathLike,
    outdir: PathLike,
    species: str = "human",
    value_threshold: float = 0,
    fasta_path: Optional[PathLike] = None,
    reference_dir: Optional[PathLike] = None,
    make_plots: bool = True,
    plot_types: Optional[Sequence[str]] = None,
    percentage_ranges: Optional[List[Dict[str, object]]] = None,
    min_percentage: int = 0,
    max_percentage: int = 100,
    sequence_logo_threshold: Optional[float] = None,
    gene_plots: Optional[Sequence[str]] = None,
    gene_plot_highlight: Optional[Union[int, tuple[int, int]]] = None,
    gene_plot_min_threshold: float = 0,
    add_all_rna_plot: bool = False,
) -> Dict[str, object]:
    """
    Run the main MitoCAPE workflow on one input file.

    Parameters
    ----------
    input_file
        Input CSV containing mutation/sample data.
    outdir
        Output directory.
    species
        One of: 'human', 'mouse', 'fly'.
    value_threshold
        Minimum threshold used in effect classification and selected plots.
    fasta_path
        Optional explicit FASTA path. If omitted, bundled/default FASTA is used.
    reference_dir
        Optional directory containing bundled reference FASTA files.
    make_plots
        Whether to generate plots.
    plot_types
        Optional subset of plots to generate. If None, a small default set is used.
    percentage_ranges
        Optional range definitions for range-based plots and summaries.
    min_percentage, max_percentage
        Used for histogram/line plots.
    sequence_logo_threshold
        Threshold for sequence-logo generation. If None, no sequence logo is made.
    gene_plots
        Optional list of genes for gene-specific plots.
    gene_plot_highlight
        Optional highlighted position or (start, end) range for gene plots.
    gene_plot_min_threshold
        Minimum percentage threshold for gene codon plots.
    add_all_rna_plot
        Whether to generate the all-RNA-type plot.

    Returns
    -------
    dict
        Dictionary with paths and summary outputs.
    """
    input_file = Path(input_file)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    ref = get_reference_data(
        species=species,
        fasta_path=fasta_path,
        reference_dir=reference_dir,
    )
    reference_df = build_reference_dataframe(
        species=species,
        fasta_path=fasta_path,
        reference_dir=reference_dir,
    )

    processed_file = process_file(
        input_file=input_file,
        reference_df=reference_df,
        codon_table=ref.codon_table,
        coding_regions=ref.coding_regions,
        output_dir=outdir,
        value_threshold=value_threshold,
    )

    outputs: Dict[str, object] = {
        "input_file": str(input_file),
        "processed_file": str(processed_file),
        "plots": [],
        "summary": None,
        "range_summary": None,
    }

    summary = calculate_effect_counts_and_percentages(processed_file)
    outputs["summary"] = summary

    if percentage_ranges:
        range_summary = calculate_effect_counts_and_percentages_ranges(
            processed_file,
            percentage_ranges=percentage_ranges,
        )
        outputs["range_summary"] = range_summary

    if not make_plots:
        return outputs

    if plot_types is None:
        plot_types = [
            "standard",
            "top_positions",
            "bubble",
        ]

    plot_paths: List[str] = []

    if "standard" in plot_types:
        plot_paths.append(str(create_standard_plot(processed_file, outdir, value_threshold=value_threshold)))

    if "top_positions" in plot_types:
        plot_paths.extend(
            str(p) for p in create_plot_most_changed_positions(
                processed_file,
                outdir,
                value_threshold=value_threshold,
            )
        )

    if "pie" in plot_types:
        plot_paths.append(str(create_pie_chart(processed_file, outdir, value_threshold=value_threshold)))

    if "donut" in plot_types:
        plot_paths.append(str(create_donut_chart(processed_file, outdir, value_threshold=value_threshold)))

    if "bubble" in plot_types:
        plot_paths.append(str(create_bubble_plot(processed_file, outdir)))

    if "histogram" in plot_types:
        plot_paths.append(
            str(
                create_histogram_per_sample(
                    processed_file,
                    outdir,
                    min_percentage=min_percentage,
                    max_percentage=max_percentage,
                )
            )
        )

    if "line" in plot_types:
        plot_paths.append(
            str(
                create_line_plot_per_sample(
                    processed_file,
                    outdir,
                    min_percentage=min_percentage,
                    max_percentage=max_percentage,
                )
            )
        )

    if add_all_rna_plot or "all_rna" in plot_types:
        plot_paths.extend(str(p) for p in create_plot_all_rna_types(processed_file, outdir))

    if sequence_logo_threshold is not None and "sequence_logo" in plot_types:
        plot_paths.append(
            str(
                create_sequence_logo_from_percentage(
                    processed_file,
                    outdir,
                    value_threshold=sequence_logo_threshold,
                )
            )
        )

    if percentage_ranges:
        if "ranges" in plot_types:
            plot_paths.append(str(create_ranges_plot(processed_file, outdir, percentage_ranges)))

        if "bubble_ranges" in plot_types:
            plot_paths.extend(
                str(p) for p in create_bubble_ranges_plot(
                    processed_file,
                    outdir,
                    percentage_ranges,
                )
            )

        if "bubble_ranges_1fig" in plot_types:
            plot_paths.append(
                str(
                    create_bubble_ranges_1fig_plot(
                        processed_file,
                        outdir,
                        percentage_ranges,
                    )
                )
            )

        if "combined_bubble" in plot_types:
            plot_paths.append(
                str(
                    create_combined_bubble_plot(
                        processed_file,
                        outdir,
                        percentage_ranges,
                    )
                )
            )

        if outputs["range_summary"]:
            range_summary_df = pd.DataFrame(outputs["range_summary"])

            if "overview_counts" in plot_types:
                if {"File", "Range", "Silent Row Count", "AA Change Row Count"}.issubset(range_summary_df.columns):
                    counts_df = range_summary_df.rename(
                        columns={
                            "Silent Row Count": "Silent Count",
                            "AA Change Row Count": "AA Change Count",
                        }
                    )
                    plot_paths.append(str(create_overview_plot_counts(counts_df, outdir, percentage_ranges)))

            if "overview_percentages" in plot_types:
                if {"File", "Range", "Silent Percentage", "AA Change Percentage"}.issubset(range_summary_df.columns):
                    plot_paths.append(
                        str(create_overview_plot_percentages(range_summary_df, outdir, percentage_ranges))
                    )

    if gene_plots:
        for gene in gene_plots:
            if "gene_mut" in plot_types:
                plot_paths.append(
                    str(
                        create_fig_gene_mut(
                            processed_file,
                            outdir,
                            gene=gene,
                            highlight=gene_plot_highlight,
                        )
                    )
                )

            if "gene_codon" in plot_types:
                plot_paths.append(
                    str(
                        create_fig_gene_codon(
                            processed_file,
                            outdir,
                            gene=gene,
                            min_threshold=gene_plot_min_threshold,
                            highlight=gene_plot_highlight,
                        )
                    )
                )

    outputs["plots"] = plot_paths
    return outputs


def enrich_with_gnomad(
    processed_file: PathLike,
    outdir: PathLike,
    gnomad_df: pd.DataFrame,
    suffix: str = "_gnomAD.csv",
) -> Dict[str, str]:
    """
    Add gnomAD annotations to an already-processed file.
    """
    processed_file = Path(processed_file)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    enriched_df = add_gnomad_metrics(
        processed_file=processed_file,
        output_dir=outdir,
        gnomad_df=gnomad_df,
        suffix=suffix,
    )

    enriched_file = outdir / f"{processed_file.stem}{suffix}"

    return {
        "processed_file": str(processed_file),
        "enriched_file": str(enriched_file),
        "rows": str(len(enriched_df)),
    }


__all__ = [
    "run_analysis",
    "enrich_with_gnomad",
]
