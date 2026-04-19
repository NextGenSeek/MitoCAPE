from __future__ import annotations

from dataclasses import dataclass
from importlib.resources import files
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Union

import pandas as pd
import Bio.Data.CodonTable
from Bio import SeqIO
from Bio.Seq import Seq


PathLike = Union[str, Path]
CodingRegion = Tuple[int, int, int, str, str]


@dataclass(frozen=True)
class ReferenceData:
    species: str
    fasta_path: Path
    gene_map_path: Path
    dna_sequence: Seq
    coding_regions: List[CodingRegion]
    codon_table: Bio.Data.CodonTable.CodonTable


SPECIES_CONFIG: Dict[str, Dict[str, object]] = {
    "human": {
        "fasta_filename": "human_chrM.fasta",
        "gene_map_filename": "human_gene_map.csv",
        "codon_table_id": 2,
    },
    "mouse": {
        "fasta_filename": "mouse_chrM.fasta",
        "gene_map_filename": "mouse_gene_map.csv",
        "codon_table_id": 2,
    },
    "fly": {
        "fasta_filename": "fly_chrM.fasta",
        "gene_map_filename": "fly_gene_map.csv",
        "codon_table_id": 5,
    },
}


def _normalise_species(species: str) -> str:
    species_norm = species.strip().lower()
    if species_norm not in SPECIES_CONFIG:
        valid = ", ".join(sorted(SPECIES_CONFIG))
        raise ValueError(f"Unsupported species '{species}'. Valid options: {valid}")
    return species_norm


def _package_reference_path(filename: str) -> Path:
    return Path(files("mitocape.data.references").joinpath(filename))


def get_default_fasta_path(species: str) -> Path:
    species_norm = _normalise_species(species)
    return _package_reference_path(str(SPECIES_CONFIG[species_norm]["fasta_filename"]))


def get_default_gene_map_path(species: str) -> Path:
    species_norm = _normalise_species(species)
    return _package_reference_path(str(SPECIES_CONFIG[species_norm]["gene_map_filename"]))


def load_fasta_sequence(fasta_path: PathLike) -> Seq:
    fasta_path = Path(fasta_path)
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    with fasta_path.open() as handle:
        record = next(SeqIO.parse(handle, "fasta"))
    return Seq(str(record.seq))


def load_gene_map(gene_map_path: PathLike) -> pd.DataFrame:
    gene_map_path = Path(gene_map_path)
    if not gene_map_path.exists():
        raise FileNotFoundError(f"Gene map file not found: {gene_map_path}")

    df = pd.read_csv(gene_map_path)
    required = {"Gene", "Start", "End", "Strand", "RNA_Type"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Gene map is missing required columns: {sorted(missing)}")

    df = df.copy()
    df["Start"] = df["Start"].astype(int)
    df["End"] = df["End"].astype(int)
    df["Strand"] = df["Strand"].astype(int)
    df["Gene"] = df["Gene"].astype(str)
    df["RNA_Type"] = df["RNA_Type"].astype(str)

    return df


def gene_map_to_coding_regions(gene_map_df: pd.DataFrame) -> List[CodingRegion]:
    return [
        (int(row.Start), int(row.End), int(row.Strand), str(row.Gene), str(row.RNA_Type))
        for row in gene_map_df.itertuples(index=False)
    ]


def get_codon_table(species: str):
    species_norm = _normalise_species(species)
    codon_table_id = int(SPECIES_CONFIG[species_norm]["codon_table_id"])
    return Bio.Data.CodonTable.unambiguous_dna_by_id[codon_table_id]


def get_reference_data(
    species: str,
    fasta_path: Optional[PathLike] = None,
    gene_map_path: Optional[PathLike] = None,
) -> ReferenceData:
    species_norm = _normalise_species(species)

    if fasta_path is None:
        fasta_path = get_default_fasta_path(species_norm)
    if gene_map_path is None:
        gene_map_path = get_default_gene_map_path(species_norm)

    fasta_path = Path(fasta_path)
    gene_map_path = Path(gene_map_path)

    dna_sequence = load_fasta_sequence(fasta_path)
    gene_map_df = load_gene_map(gene_map_path)
    coding_regions = gene_map_to_coding_regions(gene_map_df)
    codon_table = get_codon_table(species_norm)

    return ReferenceData(
        species=species_norm,
        fasta_path=fasta_path,
        gene_map_path=gene_map_path,
        dna_sequence=dna_sequence,
        coding_regions=coding_regions,
        codon_table=codon_table,
    )


def reverse_complement(sequence: str) -> str:
    return str(Seq(sequence).reverse_complement())


def translate_codon(codon: str, codon_table) -> str:
    codon = str(codon).upper()
    if len(codon) != 3 or any(base not in {"A", "C", "G", "T"} for base in codon):
        return ""

    if codon in getattr(codon_table, "stop_codons", []):
        return "X"

    return codon_table.forward_table.get(codon, "X")


def translate_with_positions(
    dna_sequence: Union[str, Seq],
    codon_table,
    regions: Sequence[CodingRegion],
) -> pd.DataFrame:
    dna_sequence = str(dna_sequence).upper()
    sequence_length = len(dna_sequence)

    original_nucleotide = list(dna_sequence)
    codons = [""] * sequence_length
    amino_acids = [""] * sequence_length
    positions_within_codon = [pd.NA] * sequence_length
    gene_names = [""] * sequence_length
    rna_types = [""] * sequence_length

    for start, end, strand, gene_name, rna_type in regions:
        if strand == 1:
            for i in range(start - 1, end, 3):
                codon = dna_sequence[i:i + 3]
                if len(codon) != 3:
                    continue

                if rna_type == "mRNA":
                    amino_acid = translate_codon(codon, codon_table)
                    for pos in range(3):
                        idx = i + pos
                        if idx < sequence_length:
                            codons[idx] = codon
                            amino_acids[idx] = amino_acid
                            positions_within_codon[idx] = pos + 1

                for pos in range(3):
                    idx = i + pos
                    if idx < sequence_length:
                        gene_names[idx] = gene_name
                        rna_types[idx] = rna_type

        elif strand == -1:
            rev_dna_sequence = reverse_complement(dna_sequence)
            for i in range(sequence_length - end, sequence_length - start + 1, 3):
                codon = rev_dna_sequence[i:i + 3]
                if len(codon) != 3:
                    continue

                if rna_type == "mRNA":
                    amino_acid = translate_codon(codon, codon_table)
                    for pos in range(3):
                        idx = sequence_length - (i + pos) - 1
                        if 0 <= idx < sequence_length:
                            codons[idx] = codon
                            amino_acids[idx] = amino_acid
                            positions_within_codon[idx] = pos + 1

                for pos in range(3):
                    idx = sequence_length - (i + pos) - 1
                    if 0 <= idx < sequence_length:
                        gene_names[idx] = gene_name
                        rna_types[idx] = rna_type
        else:
            raise ValueError(f"Invalid strand value: {strand}")

    return pd.DataFrame(
        {
            "position": range(1, sequence_length + 1),
            "Original Nucleotide": original_nucleotide,
            "CodonSequence": codons,
            "Amino Acid": amino_acids,
            "Positions within Codon": pd.array(positions_within_codon, dtype="Int64"),
            "Gene": gene_names,
            "RNA Type": rna_types,
        }
    )


def build_reference_dataframe(
    species: str,
    fasta_path: Optional[PathLike] = None,
    gene_map_path: Optional[PathLike] = None,
) -> pd.DataFrame:
    ref = get_reference_data(
        species=species,
        fasta_path=fasta_path,
        gene_map_path=gene_map_path,
    )
    return translate_with_positions(
        dna_sequence=ref.dna_sequence,
        codon_table=ref.codon_table,
        regions=ref.coding_regions,
    )


__all__ = [
    "CodingRegion",
    "ReferenceData",
    "SPECIES_CONFIG",
    "build_reference_dataframe",
    "gene_map_to_coding_regions",
    "get_codon_table",
    "get_default_fasta_path",
    "get_default_gene_map_path",
    "get_reference_data",
    "load_fasta_sequence",
    "load_gene_map",
    "reverse_complement",
    "translate_codon",
    "translate_with_positions",
]