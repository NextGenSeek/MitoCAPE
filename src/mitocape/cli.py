from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Optional

import pandas as pd

from .api import enrich_with_gnomad, run_analysis


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="mitocape",
        description="Run MitoCAPE analysis on mitochondrial mutation data.",
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Path to the input CSV file.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Directory where results will be written.",
    )
    parser.add_argument(
        "--species",
        default="human",
        choices=["human", "mouse", "fly"],
        help="Species reference to use (default: human).",
    )
    parser.add_argument(
        "--value-threshold",
        type=float,
        default=0,
        help="Minimum threshold used for effect classification and selected plots.",
    )
    parser.add_argument(
        "--fasta-path",
        default=None,
        help="Optional explicit FASTA file path. If omitted, the default bundled FASTA is used.",
    )
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Skip plot generation.",
    )
    parser.add_argument(
        "--plot-types",
        nargs="+",
        default=None,
        help=(
            "Optional list of plot types to generate. "
            "Examples: standard top_positions bubble pie donut histogram line "
            "all_rna sequence_logo ranges bubble_ranges bubble_ranges_1fig "
            "combined_bubble overview_counts overview_percentages gene_mut gene_codon"
        ),
    )

    parser.add_argument(
        "--min-percentage",
        type=int,
        default=0,
        help="Minimum percentage for histogram/line plots (default: 0).",
    )
    parser.add_argument(
        "--max-percentage",
        type=int,
        default=100,
        help="Maximum percentage for histogram/line plots (default: 100).",
    )
    parser.add_argument(
        "--sequence-logo-threshold",
        type=float,
        default=None,
        help="Threshold for sequence logo generation.",
    )

    parser.add_argument(
        "--gene-plots",
        nargs="+",
        default=None,
        help="Optional gene names for gene-specific plots.",
    )
    parser.add_argument(
        "--gene-plot-min-threshold",
        type=float,
        default=0,
        help="Minimum percentage threshold for gene codon plots.",
    )
    parser.add_argument(
        "--highlight-start",
        type=int,
        default=None,
        help="Optional start position for highlighted region in gene plots.",
    )
    parser.add_argument(
        "--highlight-end",
        type=int,
        default=None,
        help="Optional end position for highlighted region in gene plots.",
    )
    parser.add_argument(
        "--add-all-rna-plot",
        action="store_true",
        help="Generate the all-RNA-type plot.",
    )

    parser.add_argument(
        "--ranges-file",
        default=None,
        help=(
            "Optional CSV file defining percentage ranges. "
            "Expected columns: title, filename_suffix, min_percentage, max_percentage, color"
        ),
    )

    parser.add_argument(
        "--gnomad-file",
        default=None,
        help="Optional CSV file with gnomAD annotations to merge after processing.",
    )
    parser.add_argument(
        "--gnomad-suffix",
        default="_gnomAD.csv",
        help="Suffix for gnomAD-enriched output file.",
    )
    parser.add_argument(
        "--gene-map-path",
        default=None,
        help="Optional custom gene map CSV file.",
    )

    return parser.parse_args()


def build_highlight(highlight_start: Optional[int], highlight_end: Optional[int]):
    if highlight_start is None and highlight_end is None:
        return None
    if highlight_start is not None and highlight_end is None:
        return highlight_start
    if highlight_start is None and highlight_end is not None:
        return highlight_end
    return (highlight_start, highlight_end)


def load_percentage_ranges(ranges_file: Optional[str]):
    if ranges_file is None:
        return None

    ranges_path = Path(ranges_file)
    if not ranges_path.exists():
        raise FileNotFoundError(f"Ranges file not found: {ranges_path}")

    df = pd.read_csv(ranges_path)

    required = {"title", "filename_suffix", "min_percentage", "max_percentage", "color"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"Ranges file is missing required columns: {sorted(missing)}"
        )

    percentage_ranges = []
    for _, row in df.iterrows():
        min_pct = float(row["min_percentage"])
        max_pct = float(row["max_percentage"])

        percentage_ranges.append(
            {
                "title": row["title"],
                "filename_suffix": row["filename_suffix"],
                "color": row["color"],
                "condition": lambda data, lo=min_pct, hi=max_pct: (
                    (data["percentage"] >= lo) & (data["percentage"] < hi)
                ),
            }
        )

    return percentage_ranges


def main() -> None:
    args = parse_args()

    highlight = build_highlight(args.highlight_start, args.highlight_end)
    percentage_ranges = load_percentage_ranges(args.ranges_file)

    results = run_analysis(
        input_file=args.input,
        outdir=args.outdir,
        species=args.species,
        value_threshold=args.value_threshold,
        fasta_path=args.fasta_path,
        make_plots=not args.no_plots,
        plot_types=args.plot_types,
        percentage_ranges=percentage_ranges,
        min_percentage=args.min_percentage,
        max_percentage=args.max_percentage,
        sequence_logo_threshold=args.sequence_logo_threshold,
        gene_plots=args.gene_plots,
        gene_plot_highlight=highlight,
        gene_plot_min_threshold=args.gene_plot_min_threshold,
        add_all_rna_plot=args.add_all_rna_plot,
        gene_map_path=args.gene_map_path,
    )

    print("\nMitoCAPE analysis complete")
    print(f"Processed file: {results['processed_file']}")

    if results.get("summary") is not None:
        print("\nSummary:")
        for key, value in results["summary"].items():
            print(f"  {key}: {value}")

    if results.get("plots"):
        print("\nPlots:")
        for plot_path in results["plots"]:
            print(f"  {plot_path}")

    if args.gnomad_file:
        gnomad_df = pd.read_csv(args.gnomad_file)
        gnomad_results = enrich_with_gnomad(
            processed_file=results["processed_file"],
            outdir=args.outdir,
            gnomad_df=gnomad_df,
            suffix=args.gnomad_suffix,
        )

        print("\ngnomAD enrichment complete")
        print(f"Enriched file: {gnomad_results['enriched_file']}")


if __name__ == "__main__":
    main()
