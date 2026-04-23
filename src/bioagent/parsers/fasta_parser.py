"""
fasta_parser.py

Purpose: Parse FASTA files into structured Python objects.
Inputs:  Path to a FASTA file (.fasta, .fa, .fna, .faa)
Outputs: List of FastaRecord dataclass objects

Author:  Emmanuel Ogbu (Manny)
Date:    2026-04-22
"""

import logging
from pathlib import Path
from dataclasses import dataclass
from Bio import SeqIO

from bioagent.utils.logger import get_logger

logger = get_logger(__name__)


@dataclass
class FastaRecord:
    """Structured representation of a single FASTA record."""
    seq_id: str
    description: str
    sequence: str
    length: int


def parse_fasta(file_path: str | Path) -> list[FastaRecord]:
    """
    Parse a FASTA file into a list of FastaRecord objects.

    Args:
        file_path: Path to the FASTA file.

    Returns:
        List of FastaRecord objects.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If no valid records found.
    """
    path = Path(file_path)

    if not path.exists():
        raise FileNotFoundError(f"FASTA file not found: {path}")

    records = []

    try:
        # SeqIO.parse is a generator — memory efficient for large files
        for record in SeqIO.parse(path, "fasta"):
            sequence_str = str(record.seq)
            records.append(
                FastaRecord(
                    seq_id=record.id,
                    description=record.description,
                    sequence=sequence_str,
                    length=len(sequence_str),
                )
            )
    except Exception as e:
        logger.error(f"Failed to parse FASTA file {path}: {e}")
        raise

    if not records:
        raise ValueError(f"No valid FASTA records found in {path}")

    logger.info(f"Parsed {len(records)} records from {path.name}")
    return records