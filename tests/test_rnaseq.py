"""
test_rnaseq.py

Unit tests for the RNA-seq differential expression pipeline.
Author:  Emmanuel Ogbu (Manny)
Date:    2026-04-24
"""

import pytest
from pathlib import Path
from bioagent.pipelines.rnaseq import run_rnaseq_pipeline, RNAseqResult


@pytest.fixture
def counts_file() -> Path:
    return Path("data/sample/counts.csv")


def test_rnaseq_returns_result(counts_file: Path) -> None:
    result = run_rnaseq_pipeline(
        counts_file, "healthy", "cancer",
        output_dir="outputs", use_rag=False
    )
    assert isinstance(result, RNAseqResult)


def test_rnaseq_gene_counts(counts_file: Path) -> None:
    result = run_rnaseq_pipeline(
        counts_file, "healthy", "cancer",
        output_dir="outputs", use_rag=False
    )
    assert result.total_genes == 12
    assert result.genes_tested <= result.total_genes


def test_rnaseq_finds_deg(counts_file: Path) -> None:
    result = run_rnaseq_pipeline(
        counts_file, "healthy", "cancer",
        output_dir="outputs", use_rag=False
    )
    total_sig = result.upregulated_count + result.downregulated_count
    assert total_sig > 0


def test_rnaseq_generates_plots(counts_file: Path) -> None:
    result = run_rnaseq_pipeline(
        counts_file, "healthy", "cancer",
        output_dir="outputs", use_rag=False
    )
    assert len(result.plot_paths) >= 2
    for p in result.plot_paths:
        assert Path(p).exists()


def test_rnaseq_missing_file() -> None:
    with pytest.raises(FileNotFoundError):
        run_rnaseq_pipeline(
            "nonexistent.csv", "healthy", "cancer",
            use_rag=False
        )