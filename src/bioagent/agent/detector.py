"""
detector.py

Purpose: Auto-detect bioinformatics file types from uploaded files.
         This is the entry point for the agentic system — every uploaded
         file passes through here before any pipeline is selected.

Inputs:  Path to any uploaded bioinformatics file
Outputs: DetectionResult dataclass containing:
         - file_type (e.g. "FASTQ", "FASTA", "VCF")
         - confidence (0.0 to 1.0)
         - metadata (encoding, estimated read length, etc.)
         - explanation (human-readable reason for the decision)

Author:  Emmanuel Ogbu (Manny)
Date:    2026-04-22
"""

import re
import logging
from pathlib import Path
from dataclasses import dataclass, field

# Our centralised logger — imported from utils so format is consistent
from bioagent.utils.logger import get_logger

logger = get_logger(__name__)


# --- Data structures ---

@dataclass
class DetectionResult:
    """
    Structured output from the file detector.

    Using a dataclass so every downstream module gets a consistent,
    type-hinted object — not a fragile dictionary with string keys.
    """
    file_type: str              # e.g. "FASTQ", "FASTA", "VCF", "UNKNOWN"
    confidence: float           # 0.0 = no idea, 1.0 = certain
    explanation: str            # human-readable reasoning
    metadata: dict = field(default_factory=dict)  # extra info (encoding, etc.)


# --- Constants ---

# Number of lines to read from the top of the file for detection.
# We never read the whole file — some FASTQ files are 50 GB.
SAMPLE_LINES = 20


# --- Main detection function ---

def detect_file_type(file_path: str | Path) -> DetectionResult:
    """
    Detect the bioinformatics file type of an uploaded file.

    Reads only the first N lines (SAMPLE_LINES) for efficiency.
    Uses content-based detection, not just file extension — extensions
    can be wrong or missing in real-world uploads.

    Args:
        file_path: Path to the uploaded file.

    Returns:
        DetectionResult with file_type, confidence, explanation, metadata.

    Raises:
        FileNotFoundError: If the file does not exist.
        UnicodeDecodeError: If the file is binary and cannot be read as text.
    """
    path = Path(file_path)

    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")

    logger.info(f"Detecting file type for: {path.name}")

    # Read first N lines — memory safe for huge files
    try:
        sample_lines = _read_sample_lines(path, SAMPLE_LINES)
    except UnicodeDecodeError:
        # Binary file — not a standard text-based bioinformatics format
        return DetectionResult(
            file_type="UNKNOWN",
            confidence=0.0,
            explanation="File appears to be binary — cannot parse as text-based bioinformatics format."
        )

    # Run each detector in priority order — most distinctive formats first
    result = (
        _detect_fastq(sample_lines) or
        _detect_fasta(sample_lines) or
        _detect_vcf(sample_lines)   or
        _detect_gff(sample_lines)   or
        _detect_csv(sample_lines)   or
        _unknown()
    )

    logger.info(f"Detected: {result.file_type} (confidence: {result.confidence:.0%})")
    return result


# --- Private helper functions ---

def _read_sample_lines(path: Path, n: int) -> list[str]:
    """
    Read the first n non-empty lines from a file.

    Args:
        path: Path to the file.
        n: Maximum number of lines to read.

    Returns:
        List of stripped, non-empty lines.
    """
    lines = []
    # Open with UTF-8 encoding — standard for modern bioinformatics files
    with open(path, "r", encoding="utf-8", errors="strict") as f:
        for line in f:
            stripped = line.strip()
            if stripped:  # skip blank lines
                lines.append(stripped)
            if len(lines) >= n:
                break
    return lines


def _detect_fastq(lines: list[str]) -> DetectionResult | None:
    """
    Detect FASTQ format.

    FASTQ structure (repeating block of 4 lines):
      Line 1: @<sequence_id>
      Line 2: sequence (A/C/G/T/N)
      Line 3: + (separator)
      Line 4: quality scores (ASCII encoded)

    FASTQ is checked BEFORE FASTA because both start with a header line,
    but FASTQ headers start with '@' and FASTA with '>'.

    Args:
        lines: Sample lines from the file.

    Returns:
        DetectionResult if FASTQ detected, None otherwise.
    """
    if not lines:
        return None

    # FASTQ files always start with '@'
    if not lines[0].startswith("@"):
        return None

    # Check for the '+' separator line (line 3 in each block)
    has_plus = any(line == "+" or line.startswith("+") for line in lines)
    if not has_plus:
        return None

    # Try to detect Phred quality encoding
    # Line 4 (index 3) contains quality scores as ASCII characters
    encoding = "unknown"
    if len(lines) >= 4:
        quality_line = lines[3]
        encoding = _detect_phred_encoding(quality_line)

    return DetectionResult(
        file_type="FASTQ",
        confidence=0.95,
        explanation=(
            f"File starts with '@' header and contains '+' separator lines — "
            f"classic FASTQ structure. Quality encoding detected: {encoding}."
        ),
        metadata={"quality_encoding": encoding}
    )


def _detect_fasta(lines: list[str]) -> DetectionResult | None:
    """
    Detect FASTA format.

    FASTA structure:
      Line 1: ><sequence_id> [optional description]
      Line 2+: sequence (A/C/G/T/N for DNA, or amino acids for protein)

    Args:
        lines: Sample lines from the file.

    Returns:
        DetectionResult if FASTA detected, None otherwise.
    """
    if not lines:
        return None

    if not lines[0].startswith(">"):
        return None

    # Check if sequence lines look like DNA, RNA, or protein
    sequence_lines = [l for l in lines if not l.startswith(">")]
    seq_type = "unknown"
    if sequence_lines:
        seq_type = _classify_sequence(sequence_lines[0])

    return DetectionResult(
        file_type="FASTA",
        confidence=0.95,
        explanation=(
            f"File starts with '>' header line — standard FASTA format. "
            f"Sequence type appears to be: {seq_type}."
        ),
        metadata={"sequence_type": seq_type}
    )


def _detect_vcf(lines: list[str]) -> DetectionResult | None:
    """
    Detect VCF (Variant Call Format) files.

    VCF files always start with '##fileformat=VCFv' on line 1.
    Used for storing genetic variants (SNPs, indels, structural variants).

    Args:
        lines: Sample lines from the file.

    Returns:
        DetectionResult if VCF detected, None otherwise.
    """
    if not lines:
        return None

    # VCF spec requires this exact header on line 1
    if lines[0].startswith("##fileformat=VCF"):
        # Extract version number from header
        version = lines[0].split("=")[-1] if "=" in lines[0] else "unknown"
        return DetectionResult(
            file_type="VCF",
            confidence=0.99,
            explanation=f"File begins with VCF format declaration. Version: {version}.",
            metadata={"vcf_version": version}
        )

    # Some VCF files have ## meta lines without the fileformat line
    vcf_headers = [l for l in lines if l.startswith("##")]
    if len(vcf_headers) >= 3:
        return DetectionResult(
            file_type="VCF",
            confidence=0.80,
            explanation="Multiple '##' meta-information lines detected — likely VCF format.",
            metadata={}
        )

    return None


def _detect_gff(lines: list[str]) -> DetectionResult | None:
    """
    Detect GFF/GTF (General Feature Format) files.

    GFF3 files start with '##gff-version 3'.
    GTF files have tab-separated columns with specific field structure.
    Used for genome annotation — gene locations, exons, UTRs etc.

    Args:
        lines: Sample lines from the file.

    Returns:
        DetectionResult if GFF/GTF detected, None otherwise.
    """
    if not lines:
        return None

    if lines[0].startswith("##gff-version"):
        return DetectionResult(
            file_type="GFF",
            confidence=0.99,
            explanation="File begins with GFF version declaration — genome annotation format.",
            metadata={"format": "GFF3"}
        )

    # GTF detection — tab-separated with 9 columns, col 3 is feature type
    data_lines = [l for l in lines if not l.startswith("#")]
    if data_lines:
        cols = data_lines[0].split("\t")
        if len(cols) == 9:
            return DetectionResult(
                file_type="GTF",
                confidence=0.85,
                explanation="9-column tab-separated structure matches GTF genome annotation format.",
                metadata={"format": "GTF"}
            )

    return None


def _detect_csv(lines: list[str]) -> DetectionResult | None:
    """
    Detect CSV/TSV tabular data files.

    Checks for consistent delimiter (comma or tab) across sample lines.

    Args:
        lines: Sample lines from the file.

    Returns:
        DetectionResult if CSV/TSV detected, None otherwise.
    """
    if not lines:
        return None

    # Check comma-separated
    comma_counts = [line.count(",") for line in lines[:5]]
    if min(comma_counts) > 0 and max(comma_counts) == min(comma_counts):
        return DetectionResult(
            file_type="CSV",
            confidence=0.80,
            explanation=f"Consistent comma-separated columns detected ({comma_counts[0]+1} columns).",
            metadata={"delimiter": ",", "columns": comma_counts[0] + 1}
        )

    # Check tab-separated
    tab_counts = [line.count("\t") for line in lines[:5]]
    if min(tab_counts) > 0 and max(tab_counts) == min(tab_counts):
        return DetectionResult(
            file_type="TSV",
            confidence=0.80,
            explanation=f"Consistent tab-separated columns detected ({tab_counts[0]+1} columns).",
            metadata={"delimiter": "\t", "columns": tab_counts[0] + 1}
        )

    return None


def _unknown() -> DetectionResult:
    """Return a fallback result when no format is recognised."""
    return DetectionResult(
        file_type="UNKNOWN",
        confidence=0.0,
        explanation="File format could not be determined from the first 20 lines."
    )


# --- Sequence analysis helpers ---

def _classify_sequence(sequence: str) -> str:
    """
    Classify a biological sequence as DNA, RNA, or protein.

    Args:
        sequence: Raw sequence string (uppercase expected).

    Returns:
        One of: "DNA", "RNA", "protein", "unknown"
    """
    seq = sequence.upper()

    # DNA: only A, C, G, T, N (N = unknown base)
    if re.fullmatch(r"[ACGTN]+", seq):
        return "DNA"

    # RNA: same but U instead of T
    if re.fullmatch(r"[ACGUN]+", seq):
        return "RNA"

    # Protein: standard 20 amino acid single-letter codes
    if re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", seq):
        return "protein"

    return "unknown"


def _detect_phred_encoding(quality_line: str) -> str:
    """
    Detect Phred quality score encoding from a FASTQ quality line.

    Phred+33 (Sanger/Illumina 1.8+): ASCII 33-73  → chars ! to I
    Phred+64 (older Illumina):        ASCII 64-104 → chars @ to h

    This matters because misidentifying encoding causes every quality
    score to be wrong — downstream QC and trimming will be incorrect.

    Reference: Cock et al. (2010) Nucleic Acids Research
    https://doi.org/10.1093/nar/gkp1137

    Args:
        quality_line: The quality score string from a FASTQ record.

    Returns:
        "Phred+33", "Phred+64", or "unknown"
    """
    if not quality_line:
        return "unknown"

    # Get ASCII values of all quality characters
    ascii_values = [ord(c) for c in quality_line]
    min_ascii = min(ascii_values)
    max_ascii = max(ascii_values)

    # Phred+33: minimum ASCII value is 33 ('!')
    # Phred+64: minimum ASCII value is 64 ('@')
    if min_ascii < 59:
        return "Phred+33"
    elif min_ascii >= 64:
        return "Phred+64"
    else:
        return "unknown"