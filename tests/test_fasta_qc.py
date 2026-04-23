"""
test_fasta_qc.py

Unit tests for the FASTA QC pipeline.
Author:  Emmanuel Ogbu (Manny)
Date:    2026-04-23
"""

import pytest
from pathlib import Path
from bioagent.pipelines.fasta_qc import run_fasta_qc, QCResult


@pytest.fixture
def sample_fasta() -> Path:
    return Path("data/sample/test.fasta")


def test_fasta_qc_returns_result(sample_fasta: Path) -> None:
    """Pipeline should return a QCResult object."""
    result = run_fasta_qc(sample_fasta, output_dir="outputs", use_rag=False)
    assert isinstance(result, QCResult)


def test_fasta_qc_correct_sequence_count(sample_fasta: Path) -> None:
    """Should count sequences correctly."""
    result = run_fasta_qc(sample_fasta, output_dir="outputs", use_rag=False)
    assert result.total_sequences == 2


def test_fasta_qc_generates_plots(sample_fasta: Path) -> None:
    """Should generate 3 plot files."""
    result = run_fasta_qc(sample_fasta, output_dir="outputs", use_rag=False)
    assert len(result.plot_paths) == 3
    for plot_path in result.plot_paths:
        assert Path(plot_path).exists()


def test_fasta_qc_gc_content_range(sample_fasta: Path) -> None:
    """GC content should be between 0 and 100."""
    result = run_fasta_qc(sample_fasta, output_dir="outputs", use_rag=False)
    assert 0 <= result.mean_gc <= 100


def test_fasta_qc_warnings_generated(sample_fasta: Path) -> None:
    """Warnings list should never be empty."""
    result = run_fasta_qc(sample_fasta, output_dir="outputs", use_rag=False)
    assert len(result.warnings) > 0