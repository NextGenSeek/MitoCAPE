from __future__ import annotations

from dataclasses import dataclass
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
    dna_sequence: Seq
    coding_regions: List[CodingRegion]
    codon_table: Bio.Data.CodonTable.CodonTable


# ---------------------------------------------------------------------
# Species configuration
# ---------------------------------------------------------------------

SPECIES_CONFIG: Dict[str, Dict[str, object]] = {
    "human": {
        "fasta_filename": "chrM.fasta",
        "codon_table_id": 2,
        "coding_regions": [
            (576, 647, 1, "MT-TF", "tRNA"),
            (3307, 4262, 1, "MT-ND1", "mRNA"),
            (4470, 5511, 1, "MT-ND2", "mRNA"),
            (5904, 7445, 1, "MT-CO1", "mRNA"),
            (7586, 8269, 1, "MT-CO2", "mRNA"),
            (8366, 8572, 1, "MT-ATP8", "mRNA"),
            (8527, 9207, 1, "MT-ATP6", "mRNA"),
            (9207, 9990, 1, "MT-CO3", "mRNA"),
            (10059, 10404, 1, "MT-ND3", "mRNA"),
            (10470, 10766, 1, "MT-ND4L", "mRNA"),
            (10760, 12137, 1, "MT-ND4", "mRNA"),
            (12337, 14148, 1, "MT-ND5", "mRNA"),
            (14149, 14673, -1, "MT-ND6", "mRNA"),
            (14747, 15887, 1, "MT-CYB", "mRNA"),
            (576, 647, 1, "MT-TF", "tRNA"),
            (1602, 1670, 1, "MT-TV", "tRNA"),
            (1672, 1731, 1, "MT-L2", "tRNA"),
            (3230, 3304, 1, "MT-TI", "tRNA"),
            (3307, 3367, 1, "MT-TQ", "tRNA"),
            (4402, 4469, 1, "MT-TM", "tRNA"),
            (5512, 5579, 1, "MT-TW", "tRNA"),
            (7446, 7514, 1, "MT-TA", "tRNA"),
            (8261, 8291, 1, "MT-TN", "tRNA"),
            (8295, 8364, 1, "MT-TC", "tRNA"),
            (9991, 10058, 1, "MT-TG", "tRNA"),
            (10405, 10469, 1, "MT-TR", "tRNA"),
            (12138, 12206, 1, "MT-TH", "tRNA"),
            (12207, 12265, 1, "MT-TL1", "tRNA"),
            (12266, 12336, 1, "MT-TL2", "tRNA"),
            (14674, 14742, -1, "MT-TY", "tRNA"),
            (14791, 14858, 1, "MT-TS1", "tRNA"),
            (15888, 15953, 1, "MT-TT", "tRNA"),
            (15956, 16023, 1, "MT-TP", "tRNA"),
            (14180, 14238, -1, "MT-TG", "tRNA"),
            (14258, 14286, 1, "MT-TL", "tRNA"),
            (648, 1601, 1, "MT-RNR1", "rRNA"),
            (1671, 3229, 1, "MT-RNR2", "rRNA"),
        ],
    },
    "mouse": {
        "fasta_filename": "Mus_musculus.GRCm39.dna.chromosome.MT.fa",
        "codon_table_id": 2,
        "coding_regions": [
            (2751, 3707, 1, "MT-ND1", "mRNA"),
            (3914, 4951, 1, "MT-ND2", "mRNA"),
            (5328, 6872, 1, "MT-CO1", "mRNA"),
            (7013, 7696, 1, "MT-CO2", "mRNA"),
            (7766, 7969, 1, "MT-ATP8", "mRNA"),
            (7927, 8607, 1, "MT-ATP6", "mRNA"),
            (8607, 9390, 1, "MT-CO3", "mRNA"),
            (9459, 9806, 1, "MT-ND3", "mRNA"),
            (9877, 10173, 1, "MT-ND4L", "mRNA"),
            (10167, 11544, 1, "MT-ND4", "mRNA"),
            (11742, 13565, 1, "MT-ND5", "mRNA"),
            (13552, 14070, -1, "MT-ND6", "mRNA"),
            (14145, 15288, 1, "MT-CYB", "mRNA"),
            (68, 136, 1, "MT-TF", "tRNA"),
            (1445, 1517, 1, "MT-TV", "tRNA"),
            (1602, 1670, 1, "MT-TM", "tRNA"),
            (3230, 3304, 1, "MT-TI", "tRNA"),
            (3307, 3367, 1, "MT-TQ", "tRNA"),
            (5512, 5579, 1, "MT-TW", "tRNA"),
            (7446, 7514, 1, "MT-TA", "tRNA"),
            (8305, 8374, 1, "MT-TN", "tRNA"),
            (8610, 8676, 1, "MT-TC", "tRNA"),
            (10405, 10469, 1, "MT-TR", "tRNA"),
            (12138, 12206, 1, "MT-TH", "tRNA"),
            (12207, 12265, 1, "MT-TL1", "tRNA"),
            (12266, 12336, 1, "MT-TL2", "tRNA"),
            (14674, 14742, -1, "MT-TY", "tRNA"),
            (14791, 14858, 1, "MT-TS1", "tRNA"),
            (15888, 15953, 1, "MT-TT", "tRNA"),
            (15956, 16023, 1, "MT-TP", "tRNA"),
            (1410, 2207, 1, "MT-RNR1", "rRNA"),
            (1671, 3229, 1, "MT-RNR2", "rRNA"),
        ],
    },
    "fly": {
        "fasta_filename": "DrosMchrM.fasta",
        "codon_table_id": 5,
        "coding_regions": [
            (1, 65, 1, "MT-TI", "tRNA"),
            (97, 165, -1, "MT-TG", "tRNA"),
            (171, 239, 1, "MT-TM", "tRNA"),
            (240, 1263, 1, "MT-ND2", "mRNA"),
            (1264, 1329, 1, "MT-TT", "tRNA"),
            (1322, 1383, -1, "MT-TC", "tRNA"),
            (1403, 1468, -1, "MT-TT", "tRNA"),
            (1474, 3009, 1, "MT-COX1", "mRNA"),
            (3012, 3077, 1, "MT-TL", "tRNA"),
            (3083, 3767, 1, "MT-COX2", "mRNA"),
            (3768, 3838, 1, "MT-TL", "tRNA"),
            (3840, 3906, 1, "MT-TA", "tRNA"),
            (3907, 4068, 1, "MT-ATP8", "mRNA"),
            (4062, 4736, 1, "MT-ATP6", "mRNA"),
            (4736, 5524, 1, "MT-COX3", "mRNA"),
            (5543, 5607, 1, "MT-TG", "tRNA"),
            (5608, 5961, 1, "MT-ND3", "mRNA"),
            (5983, 6047, 1, "MT-TA", "tRNA"),
            (6062, 6125, 1, "MT-TA", "tRNA"),
            (6126, 6190, 1, "MT-TA", "tRNA"),
            (6191, 6258, 1, "MT-TS", "tRNA"),
            (6259, 6325, 1, "MT-TG", "tRNA"),
            (6344, 6408, -1, "MT-TP", "tRNA"),
            (6409, 8125, -1, "MT-ND5", "mRNA"),
            (8141, 8206, -1, "MT-TH", "tRNA"),
            (8207, 9545, -1, "MT-ND4", "mRNA"),
            (9545, 9835, -1, "MT-ND4L", "mRNA"),
            (9838, 9903, 1, "MT-TT", "tRNA"),
            (9904, 9964, -1, "MT-TP", "tRNA"),
            (9971, 10495, 1, "MT-ND6", "mRNA"),
            (10499, 11635, 1, "MT-CYTB", "mRNA"),
            (11638, 11703, 1, "MT-TS", "tRNA"),
            (11721, 12659, -1, "MT-ND1", "mRNA"),
            (12670, 12734, -1, "MT-TL", "tRNA"),
            (12735, 14058, -1, "MT-RNR", "rRNA"),
            (14058, 14130, -1, "MT-TV", "tRNA"),
            (14131, 14916, -1, "MT-RNR", "rRNA"),
        ],
    },
}


# ---------------------------------------------------------------------
# Core reference loading
# ---------------------------------------------------------------------

def _normalise_species(species: str) -> str:
    species_norm = species.strip().lower()
    if species_norm not in SPECIES_CONFIG:
        valid = ", ".join(sorted(SPECIES_CONFIG))
        raise ValueError(f"Unsupported species '{species}'. Valid options: {valid}")
    return species_norm


def load_fasta_sequence(fasta_path: PathLike) -> Seq:
    """
    Load the first sequence from a FASTA file and return it as a Bio.Seq.Seq object.
    """
    fasta_path = Path(fasta_path)
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    with fasta_path.open() as handle:
        record = next(SeqIO.parse(handle, "fasta"))
    return Seq(str(record.seq))


def get_coding_regions(species: str) -> List[CodingRegion]:
    """
    Return the coding regions for a supported species.
    """
    species_norm = _normalise_species(species)
    return list(SPECIES_CONFIG[species_norm]["coding_regions"])


def get_codon_table(species: str):
    """
    Return the BioPython codon table for a supported species.
    """
    species_norm = _normalise_species(species)
    codon_table_id = int(SPECIES_CONFIG[species_norm]["codon_table_id"])
    return Bio.Data.CodonTable.unambiguous_dna_by_id[codon_table_id]


def get_default_fasta_path(species: str, reference_dir: Optional[PathLike] = None) -> Path:
    """
    Return the default FASTA path for a species.

    If reference_dir is not provided, defaults to:
    <repo>/data/references relative to this module.
    """
    species_norm = _normalise_species(species)

    if reference_dir is None:
        # src/mitocape/references.py -> repo root / data / references
        reference_dir = Path(__file__).resolve().parents[2] / "data" / "references"

    reference_dir = Path(reference_dir)
    fasta_filename = str(SPECIES_CONFIG[species_norm]["fasta_filename"])
    return reference_dir / fasta_filename


def get_reference_data(
    species: str,
    fasta_path: Optional[PathLike] = None,
    reference_dir: Optional[PathLike] = None,
) -> ReferenceData:
    """
    Load all core reference data for a species.

    Parameters
    ----------
    species
        One of: 'human', 'mouse', 'fly'
    fasta_path
        Optional explicit FASTA path. If omitted, the default species FASTA is used.
    reference_dir
        Optional directory containing bundled reference FASTA files.

    Returns
    -------
    ReferenceData
    """
    species_norm = _normalise_species(species)

    if fasta_path is None:
        fasta_path = get_default_fasta_path(species_norm, reference_dir=reference_dir)

    fasta_path = Path(fasta_path)
    dna_sequence = load_fasta_sequence(fasta_path)
    coding_regions = get_coding_regions(species_norm)
    codon_table = get_codon_table(species_norm)

    return ReferenceData(
        species=species_norm,
        fasta_path=fasta_path,
        dna_sequence=dna_sequence,
        coding_regions=coding_regions,
        codon_table=codon_table,
    )


# ---------------------------------------------------------------------
# Sequence utilities
# ---------------------------------------------------------------------

def reverse_complement(sequence: str) -> str:
    """
    Return the reverse complement of a DNA sequence.
    """
    return str(Seq(sequence).reverse_complement())


def translate_codon(codon: str, codon_table) -> str:
    """
    Translate a codon using the supplied codon table.

    Returns:
    - amino acid single-letter code for standard codons
    - 'X' for stop codons or unknown/invalid codons
    """
    codon = str(codon).upper()
    if len(codon) != 3 or any(base not in {"A", "C", "G", "T"} for base in codon):
        return ""

    if codon in getattr(codon_table, "stop_codons", []):
        return "X"

    return codon_table.forward_table.get(codon, "X")


# ---------------------------------------------------------------------
# Reference annotation table generation
# ---------------------------------------------------------------------

def translate_with_positions(
    dna_sequence: Union[str, Seq],
    codon_table,
    regions: Sequence[CodingRegion],
) -> pd.DataFrame:
    """
    Build a nucleotide-level reference annotation table.

    Output columns:
    - position (1-based)
    - Original Nucleotide
    - CodonSequence
    - Amino Acid
    - Positions within Codon
    - Gene
    - RNA Type

    Notes
    -----
    Only mRNA regions receive codon/amino-acid/position-within-codon annotation.
    All region types receive Gene and RNA Type annotation.
    """
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
            raise ValueError(
                f"Invalid strand value {strand} for region "
                f"({start}, {end}, {strand}, {gene_name}, {rna_type})"
            )

    df = pd.DataFrame(
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

    return df


def build_reference_dataframe(
    species: str,
    fasta_path: Optional[PathLike] = None,
    reference_dir: Optional[PathLike] = None,
) -> pd.DataFrame:
    """
    Convenience wrapper: load the reference data for a species and return the
    annotated nucleotide-level reference DataFrame.
    """
    ref = get_reference_data(
        species=species,
        fasta_path=fasta_path,
        reference_dir=reference_dir,
    )
    return translate_with_positions(
        dna_sequence=ref.dna_sequence,
        codon_table=ref.codon_table,
        regions=ref.coding_regions,
    )


# ---------------------------------------------------------------------
# Optional SIFT support
# ---------------------------------------------------------------------

def load_sift_file(sift_path: PathLike, **read_csv_kwargs) -> pd.DataFrame:
    """
    Load a SIFT annotation file.

    This is intentionally generic because different SIFT exports can have
    slightly different delimiters/columns. By default it assumes tab-separated.
    """
    sift_path = Path(sift_path)
    if not sift_path.exists():
        raise FileNotFoundError(f"SIFT file not found: {sift_path}")

    default_kwargs = {"sep": "\t"}
    default_kwargs.update(read_csv_kwargs)
    return pd.read_csv(sift_path, **default_kwargs)


__all__ = [
    "CodingRegion",
    "ReferenceData",
    "SPECIES_CONFIG",
    "build_reference_dataframe",
    "get_coding_regions",
    "get_codon_table",
    "get_default_fasta_path",
    "get_reference_data",
    "load_fasta_sequence",
    "load_sift_file",
    "reverse_complement",
    "translate_codon",
    "translate_with_positions",
]
