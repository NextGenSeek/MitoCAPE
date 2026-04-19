from __future__ import annotations

from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Union

import pandas as pd

from .references import reverse_complement, translate_codon


PathLike = Union[str, Path]


# ---------------------------------------------------------------------
# Core annotation helpers
# ---------------------------------------------------------------------

def get_strand_orientation(position: int, coding_regions: Sequence[tuple]) -> Optional[int]:
    """
    Determine strand orientation for a genomic position.

    Expects coding-region tuples in the format:
    (start, end, strand, gene, rna_type)
    """
    for region in coding_regions:
        start, end, strand = region[0], region[1], region[2]
        if start <= position <= end:
            return strand
    return None


def construct_new_codon_sequence(row: pd.Series, combined_df: pd.DataFrame) -> str:
    """
    Construct the mutated codon sequence for a given row based on:
    - position within codon
    - original nucleotide
    - mutated nucleotide
    - strand orientation

    Returns an empty string if a complete codon cannot be reconstructed.
    """
    pos_in_codon = row.get("Positions within Codon")
    original_nt = row.get("Original Nucleotide", "")
    new_nt = str(row.get("NewSequence", ""))
    strand_orientation = row.get("Strand Orientation")

    if pd.isna(pos_in_codon) or strand_orientation not in (1, -1):
        return ""

    try:
        pos_in_codon = int(pos_in_codon)
    except Exception:
        return ""

    def get_nucleotide(target_position: int, codon_pos: int) -> Optional[str]:
        subset = combined_df[
            (combined_df["position"] == target_position)
            & (combined_df["Positions within Codon"] == codon_pos)
        ]
        if subset.empty:
            return None
        return subset.iloc[0]["Original Nucleotide"]

    new_codon_seq = ""

    if strand_orientation == 1:
        if pos_in_codon == 1:
            nt2 = get_nucleotide(row["position"] + 1, 2)
            nt3 = get_nucleotide(row["position"] + 2, 3)
            if nt2 is not None and nt3 is not None:
                new_codon_seq = f"{new_nt}{nt2}{nt3}" if new_nt != original_nt else f"{original_nt}{nt2}{nt3}"

        elif pos_in_codon == 2:
            nt1 = get_nucleotide(row["position"] - 1, 1)
            nt3 = get_nucleotide(row["position"] + 1, 3)
            if nt1 is not None and nt3 is not None:
                new_codon_seq = f"{nt1}{new_nt}{nt3}" if new_nt != original_nt else f"{nt1}{original_nt}{nt3}"

        elif pos_in_codon == 3:
            nt1 = get_nucleotide(row["position"] - 2, 1)
            nt2 = get_nucleotide(row["position"] - 1, 2)
            if nt1 is not None and nt2 is not None:
                new_codon_seq = f"{nt1}{nt2}{new_nt}" if new_nt != original_nt else f"{nt1}{nt2}{original_nt}"

    elif strand_orientation == -1:
        if pos_in_codon == 1:
            nt2 = get_nucleotide(row["position"] - 1, 2)
            nt3 = get_nucleotide(row["position"] - 2, 3)
            if nt2 is not None and nt3 is not None:
                new_codon_seq = f"{nt3}{nt2}{new_nt}" if new_nt != original_nt else f"{nt3}{nt2}{original_nt}"

        elif pos_in_codon == 2:
            nt1 = get_nucleotide(row["position"] + 1, 1)
            nt3 = get_nucleotide(row["position"] - 1, 3)
            if nt1 is not None and nt3 is not None:
                new_codon_seq = f"{nt3}{new_nt}{nt1}" if new_nt != original_nt else f"{nt3}{original_nt}{nt1}"

        elif pos_in_codon == 3:
            nt1 = get_nucleotide(row["position"] + 2, 1)
            nt2 = get_nucleotide(row["position"] + 1, 2)
            if nt1 is not None and nt2 is not None:
                new_codon_seq = f"{new_nt}{nt2}{nt1}" if new_nt != original_nt else f"{original_nt}{nt2}{nt1}"

        if new_codon_seq:
            new_codon_seq = reverse_complement(new_codon_seq)

    return new_codon_seq


def add_protein_positions(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add a 'PositionProtein' column indicating codon number within each mRNA gene.

    Each set of three nucleotides belonging to the same codon shares the same
    PositionProtein value.
    """
    required = {"Positions within Codon", "Gene", "RNA Type", "position"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns for add_protein_positions(): {sorted(missing)}")

    df = df.copy()
    df["PositionProtein"] = pd.NA

    current_protein_positions: Dict[str, int] = {}

    df = df.sort_values(by=["Gene", "position"]).reset_index(drop=True)

    for index, row in df.iterrows():
        if row["RNA Type"] != "mRNA" or pd.isna(row["Positions within Codon"]):
            continue

        gene = row["Gene"]
        pos_in_codon = row["Positions within Codon"]

        if pd.isna(gene) or gene == "":
            continue

        if int(pos_in_codon) == 1:
            if gene not in current_protein_positions:
                current_protein_positions[gene] = 1
            else:
                current_protein_positions[gene] += 1

        if gene in current_protein_positions:
            df.at[index, "PositionProtein"] = current_protein_positions[gene]
        else:
            current_protein_positions[gene] = 1
            df.at[index, "PositionProtein"] = 1

    df["PositionProtein"] = pd.array(df["PositionProtein"], dtype="Int64")
    return df


def determine_effect(row: pd.Series, value_threshold: Optional[float] = None) -> str:
    """
    Determine mutation effect category for a row.

    Returns:
    - 'silent'
    - 'AA_Change'
    - '' (empty string) if below threshold or not applicable
    """
    threshold = row.get("value_threshold", value_threshold)
    if threshold is None:
        raise ValueError("Provide a scalar `value_threshold` or a 'value_threshold' column.")

    percentage = pd.to_numeric(row.get("percentage"), errors="coerce")
    if pd.isna(percentage) or percentage <= threshold:
        return ""

    if row.get("Original Nucleotide") != row.get("NewSequence"):
        aa = row.get("Amino Acid", "")
        new_aa = row.get("New_AminoAcid", "")

        if aa == new_aa and aa != "":
            return "silent"
        if aa != new_aa:
            return "AA_Change"

    return ""


# ---------------------------------------------------------------------
# Main processing
# ---------------------------------------------------------------------

def process_dataframe(
    sample_df: pd.DataFrame,
    reference_df: pd.DataFrame,
    codon_table,
    coding_regions: Sequence[tuple],
    value_threshold: float = 0,
) -> pd.DataFrame:
    """
    Process a mutation/sample DataFrame by merging with reference annotation and
    calculating codon-level and amino-acid-level effects.

    Expected sample input columns:
    - position
    - ref_base
    - MutBase
    - percentage (recommended for determine_effect)

    Returns a processed DataFrame.
    """
    required = {"position", "ref_base", "MutBase"}
    missing = required - set(sample_df.columns)
    if missing:
        raise ValueError(f"Missing required columns in sample_df: {sorted(missing)}")

    df = sample_df.copy()

    # Fill missing mutated sequence with reference base
    df["NewSequence"] = df["MutBase"].fillna(df["ref_base"])

    combined_df = pd.merge(df, reference_df, on="position", how="left")

    if "Positions within Codon" not in combined_df.columns:
        raise ValueError(
            "Column 'Positions within Codon' missing after merge with reference annotation."
        )

    combined_df["Strand Orientation"] = combined_df["position"].apply(
        lambda pos: get_strand_orientation(pos, coding_regions)
    )

    combined_df["NewCodonSequence"] = combined_df.apply(
        lambda row: construct_new_codon_sequence(row, combined_df),
        axis=1,
    )

    combined_df["New_AminoAcid"] = combined_df.apply(
        lambda row: translate_codon(row["NewCodonSequence"], codon_table)
        if row.get("Amino Acid", "") != ""
        else "",
        axis=1,
    )

    if {"Positions within Codon", "Gene", "RNA Type"}.issubset(combined_df.columns):
        combined_df = add_protein_positions(combined_df)

    combined_df["effect"] = combined_df.apply(
        determine_effect,
        axis=1,
        value_threshold=value_threshold,
    )

    # Standardise stop codons as X
    combined_df.loc[
        combined_df["NewCodonSequence"].isin(["TAA", "TAG", "AGA", "AGG"]),
        "New_AminoAcid",
    ] = "X"

    return combined_df


def process_file(
    input_file: PathLike,
    reference_df: pd.DataFrame,
    codon_table,
    coding_regions: Sequence[tuple],
    output_dir: PathLike,
    value_threshold: float = 0,
    output_suffix: str = "_processed.csv",
) -> Path:
    """
    Process a mutation file and save the processed output CSV.

    Returns the output file path.
    """
    input_file = Path(input_file)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    sample_df = pd.read_csv(input_file)

    processed_df = process_dataframe(
        sample_df=sample_df,
        reference_df=reference_df,
        codon_table=codon_table,
        coding_regions=coding_regions,
        value_threshold=value_threshold,
    )

    output_file = output_dir / f"{input_file.stem}{output_suffix}"
    processed_df.to_csv(output_file, index=False)

    return output_file


# ---------------------------------------------------------------------
# Summaries
# ---------------------------------------------------------------------

def calculate_effect_counts_and_percentages(filename: PathLike) -> Optional[Dict[str, float]]:
    """
    Calculate silent vs amino-acid-changing mutation counts and percentages from a
    processed file.
    """
    filename = Path(filename)
    combined_df = pd.read_csv(filename)

    required_columns = [
        "Positions within Codon",
        "Amino Acid",
        "New_AminoAcid",
        "percentage",
        "effect",
    ]
    if not all(col in combined_df.columns for col in required_columns):
        return None

    filtered_df = combined_df.dropna(subset=["Positions within Codon", "percentage", "Amino Acid"])

    if filtered_df.empty:
        return None

    silent_count = int((filtered_df["effect"] == "silent").sum())
    aa_change_count = int((filtered_df["effect"] == "AA_Change").sum())
    total_count = int(len(filtered_df))

    silent_percentage = (silent_count / total_count) * 100 if total_count > 0 else 0.0
    aa_change_percentage = (aa_change_count / total_count) * 100 if total_count > 0 else 0.0

    return {
        "File": filename.name.replace("_processed.csv", ""),
        "Silent Count": silent_count,
        "AA Change Count": aa_change_count,
        "Silent Percentage": silent_percentage,
        "AA Change Percentage": aa_change_percentage,
    }


def calculate_effect_counts_and_percentages_ranges(
    filename: PathLike,
    percentage_ranges: List[Dict[str, object]],
) -> List[Dict[str, object]]:
    """
    Calculate range-specific effect summaries from a processed file.

    Each percentage_ranges item should contain:
    - 'title': display name
    - 'condition': callable taking a DataFrame and returning a boolean mask
    """
    filename = Path(filename)
    df = pd.read_csv(filename)

    required = [
        "Positions within Codon",
        "Amino Acid",
        "New_AminoAcid",
        "percentage",
        "effect",
        "TotalNucleotideCoverage",
    ]
    if not all(col in df.columns for col in required):
        return []

    filtered = df.dropna(subset=["Positions within Codon", "percentage", "Amino Acid"])
    if filtered.empty:
        return []

    sample_col = "Sample" if "Sample" in filtered.columns else "_SYNTH_SAMPLE_"
    if sample_col == "_SYNTH_SAMPLE_":
        filtered = filtered.copy()
        filtered[sample_col] = filename.name.replace("_processed.csv", "")

    filtered["_cov_"] = pd.to_numeric(filtered["TotalNucleotideCoverage"], errors="coerce").fillna(0)

    results: List[Dict[str, object]] = []

    for sample_id, df_sample in filtered.groupby(sample_col, dropna=False):
        total_reads = float(df_sample["_cov_"].sum())

        for rng in percentage_ranges:
            condition = rng.get("condition")
            if not callable(condition):
                continue

            try:
                raw_mask = condition(df_sample)
            except Exception:
                raw_mask = False

            mask = pd.Series(raw_mask, index=df_sample.index)
            mask = mask.reindex(df_sample.index).fillna(False).astype(bool)

            df_range = df_sample.loc[mask]

            silent_rows = int((df_range["effect"] == "silent").sum())
            aa_rows = int((df_range["effect"] == "AA_Change").sum())

            silent_reads = float(df_range.loc[df_range["effect"] == "silent", "_cov_"].sum())
            aa_reads = float(df_range.loc[df_range["effect"] == "AA_Change", "_cov_"].sum())

            denom = total_reads
            silent_pct = (silent_reads / denom * 100.0) if denom > 0 else 0.0
            aa_pct = (aa_reads / denom * 100.0) if denom > 0 else 0.0

            results.append(
                {
                    "File": filename.name.replace("_processed.csv", ""),
                    "Sample": sample_id,
                    "Range": rng.get("title", ""),
                    "All Reads Total Count": int(denom) if float(denom).is_integer() else denom,
                    "Silent Row Count": silent_rows,
                    "AA Change Row Count": aa_rows,
                    "Silent Percentage": silent_pct,
                    "AA Change Percentage": aa_pct,
                }
            )

    return results


# ---------------------------------------------------------------------
# Optional enrichment
# ---------------------------------------------------------------------

def add_gnomad_metrics(
    processed_file: PathLike,
    output_dir: PathLike,
    gnomad_df: pd.DataFrame,
    suffix: str = "_gnomAD.csv",
) -> pd.DataFrame:
    """
    Merge gnomAD-style annotations into a processed file and save the enriched output.

    Required gnomAD columns:
    - position
    - ref_base
    - MutBase
    - AF_hom
    - AF_het
    - max_observed_heteroplasmy
    """
    processed_file = Path(processed_file)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    proc = pd.read_csv(processed_file)

    required = {"position", "ref_base", "MutBase", "AF_hom", "AF_het", "max_observed_heteroplasmy"}
    missing = required - set(gnomad_df.columns)
    if missing:
        raise ValueError(f"Missing required columns in gnomad_df: {sorted(missing)}")

    enriched = proc.merge(
        gnomad_df[
            [
                "position",
                "ref_base",
                "MutBase",
                "AF_hom",
                "AF_het",
                "max_observed_heteroplasmy",
            ]
        ],
        on=["position", "ref_base", "MutBase"],
        how="left",
    )

    output_file = output_dir / f"{processed_file.stem}{suffix}"
    enriched.to_csv(output_file, index=False)

    return enriched


__all__ = [
    "add_gnomad_metrics",
    "add_protein_positions",
    "calculate_effect_counts_and_percentages",
    "calculate_effect_counts_and_percentages_ranges",
    "construct_new_codon_sequence",
    "determine_effect",
    "get_strand_orientation",
    "process_dataframe",
    "process_file",
]
