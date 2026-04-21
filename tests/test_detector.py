"""
test_detector.py

Purpose: Unit tests for the file type auto-detection module.
         Tests each format detector independently and together.

Author:  Emmanuel Ogbu (Manny)
Date:    2026-04-22
"""

import pytest
from pathlib import Path
from bioagent.agent.detector import detect_file_type, DetectionResult


# --- Fixtures ---

@pytest.fixture
def sample_dir() -> Path:
    """Path to the sample data directory."""
    return Path(__file__).parent.parent / "data" / "sample"

@pytest.fixture
def fasta_file(sample_dir) -> Path:
    return sample_dir / "test.fasta"

@pytest.fixture
def fastq_file(sample_dir) -> Path:
    return sample_dir / "test.fastq"


# --- Tests ---

def test_detect_fasta(fasta_file: Path) -> None:
    result = detect_file_type(fasta_file)
    assert result.file_type == "FASTA"
    assert result.confidence >= 0.9

def test_detect_fastq(fastq_file: Path) -> None:
    result = detect_file_type(fastq_file)
    assert result.file_type == "FASTQ"
    assert result.confidence >= 0.9

def test_detect_missing_file() -> None:
    with pytest.raises(FileNotFoundError):
        detect_file_type("nonexistent.fasta")

def test_result_is_detection_result(fasta_file: Path) -> None:
    result = detect_file_type(fasta_file)
    assert isinstance(result, DetectionResult)
    assert result.explanation != ""