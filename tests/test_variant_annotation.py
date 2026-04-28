"""
test_variant_annotation.py

Unit tests for the variant annotation pipeline.
Author:  Emmanuel Ogbu (Manny)
Date:    2026-04-24
"""

import pytest
from pathlib import Path
from bioagent.pipelines.variant_annotation import (
    run_variant_pipeline, VariantResult
)


@pytest.fixture
def vcf_file() -> Path:
    return Path("data/sample/test.vcf")


def test_variant_pipeline_returns_result(vcf_file: Path) -> None:
    result = run_variant_pipeline(vcf_file, output_dir="outputs", use_rag=False)
    assert isinstance(result, VariantResult)


def test_variant_pipeline_parses_all(vcf_file: Path) -> None:
    result = run_variant_pipeline(vcf_file, output_dir="outputs", use_rag=False)
    assert result.total_variants == 8


def test_variant_pipeline_finds_pathogenic(vcf_file: Path) -> None:
    result = run_variant_pipeline(vcf_file, output_dir="outputs", use_rag=False)
    assert result.pathogenic_count > 0


def test_variant_pipeline_generates_plots(vcf_file: Path) -> None:
    result = run_variant_pipeline(vcf_file, output_dir="outputs", use_rag=False)
    assert len(result.plot_paths) >= 2
    for p in result.plot_paths:
        assert Path(p).exists()


def test_variant_pipeline_missing_file() -> None:
    with pytest.raises(FileNotFoundError):
        run_variant_pipeline("nonexistent.vcf", use_rag=False)