"""
router.py

Purpose: Route uploaded files to the correct bioinformatics pipeline.
         This is the decision-making module of the agentic system.
         It receives a DetectionResult from detector.py and returns
         the appropriate pipeline function to run.

         Decision logic:
         - FASTA  → FASTA QC pipeline
         - FASTQ  → FASTA QC pipeline (quality metrics still apply)
         - VCF    → Variant Annotation pipeline
         - CSV    → RNA-seq pipeline (assumes count matrix)
         - TSV    → RNA-seq pipeline (assumes count matrix)
         - UNKNOWN → raises error with helpful message

Inputs:  DetectionResult from detector.py + file path
Outputs: PipelineResult (one of QCResult, RNAseqResult, VariantResult)

Author:  Emmanuel Ogbu (Manny)
Date:    2026-04-28
"""

from pathlib import Path
from dataclasses import dataclass
from typing import Any

from bioagent.agent.detector import detect_file_type, DetectionResult
from bioagent.utils.logger import get_logger

logger = get_logger(__name__)


@dataclass
class RoutingDecision:
    """
    Records what the router decided and why.
    Stored alongside results so the frontend can explain the decision.
    """
    file_type: str          # detected file type
    pipeline_name: str      # which pipeline was selected
    confidence: float       # detection confidence (0-1)
    reasoning: str          # human-readable explanation of the decision


def route_file(
    file_path: str | Path,
    output_dir: str | Path = "outputs",
    use_rag: bool = True,
    rnaseq_control: str = "control",
    rnaseq_treatment: str = "treatment"
) -> tuple[Any, RoutingDecision]:
    """
    Auto-detect a file's type and route it to the correct pipeline.

    This is the main entry point for the agentic system. A user uploads
    any bioinformatics file and this function handles everything:
    detection → routing → pipeline execution → results.

    Args:
        file_path: Path to the uploaded file.
        output_dir: Directory to save output plots.
        use_rag: Whether to use RAG for biological interpretation.
        rnaseq_control: Column prefix for control samples in CSV files.
        rnaseq_treatment: Column prefix for treatment samples in CSV files.

    Returns:
        Tuple of (pipeline_result, RoutingDecision).
        pipeline_result is one of: QCResult, RNAseqResult, VariantResult.

    Raises:
        FileNotFoundError: If the file doesn't exist.
        ValueError: If the file type cannot be handled.
    """
    path = Path(file_path)

    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")

    logger.info(f"Routing file: {path.name}")

    # Step 1 — Detect file type
    detection = detect_file_type(path)
    logger.info(
        f"Detected: {detection.file_type} "
        f"(confidence: {detection.confidence:.0%})"
    )

    # Step 2 — Route to correct pipeline
    file_type = detection.file_type.upper()

    if file_type in ("FASTA", "FASTQ"):
        return _run_fasta_pipeline(path, detection, output_dir, use_rag)

    elif file_type == "VCF":
        return _run_variant_pipeline(path, detection, output_dir, use_rag)

    elif file_type in ("CSV", "TSV"):
        return _run_rnaseq_pipeline(
            path, detection, output_dir, use_rag,
            rnaseq_control, rnaseq_treatment
        )

    elif file_type == "UNKNOWN":
        raise ValueError(
            f"Could not determine file type for {path.name}. "
            f"Supported formats: FASTA, FASTQ, VCF, CSV, TSV. "
            f"Please check your file format."
        )

    else:
        raise ValueError(
            f"No pipeline available for file type: {file_type}. "
            f"Supported: FASTA, FASTQ, VCF, CSV, TSV."
        )


def _run_fasta_pipeline(
    path: Path,
    detection: DetectionResult,
    output_dir: str | Path,
    use_rag: bool
) -> tuple[Any, RoutingDecision]:
    """Route FASTA/FASTQ files to the FASTA QC pipeline."""

    # Import here to avoid circular imports and keep startup fast
    from bioagent.pipelines.fasta_qc import run_fasta_qc

    reasoning = (
        f"File detected as {detection.file_type} with "
        f"{detection.confidence:.0%} confidence. "
        f"{detection.explanation} "
        f"Routing to FASTA QC pipeline for sequence quality analysis: "
        f"GC content, length distribution, complexity assessment."
    )

    decision = RoutingDecision(
        file_type=detection.file_type,
        pipeline_name="FASTA QC Pipeline",
        confidence=detection.confidence,
        reasoning=reasoning
    )

    logger.info(f"Routing to FASTA QC pipeline.")
    result = run_fasta_qc(path, output_dir=output_dir, use_rag=use_rag)
    return result, decision


def _run_variant_pipeline(
    path: Path,
    detection: DetectionResult,
    output_dir: str | Path,
    use_rag: bool
) -> tuple[Any, RoutingDecision]:
    """Route VCF files to the variant annotation pipeline."""

    from bioagent.pipelines.variant_annotation import run_variant_pipeline

    reasoning = (
        f"File detected as VCF (Variant Call Format) with "
        f"{detection.confidence:.0%} confidence. "
        f"{detection.explanation} "
        f"Routing to Variant Annotation pipeline: quality filtering, "
        f"gene annotation, clinical significance assessment."
    )

    decision = RoutingDecision(
        file_type=detection.file_type,
        pipeline_name="Variant Annotation Pipeline",
        confidence=detection.confidence,
        reasoning=reasoning
    )

    logger.info(f"Routing to Variant Annotation pipeline.")
    result = run_variant_pipeline(path, output_dir=output_dir, use_rag=use_rag)
    return result, decision


def _run_rnaseq_pipeline(
    path: Path,
    detection: DetectionResult,
    output_dir: str | Path,
    use_rag: bool,
    control_label: str,
    treatment_label: str
) -> tuple[Any, RoutingDecision]:
    """Route CSV/TSV files to the RNA-seq pipeline."""

    from bioagent.pipelines.rnaseq import run_rnaseq_pipeline

    reasoning = (
        f"File detected as {detection.file_type} tabular data with "
        f"{detection.confidence:.0%} confidence. "
        f"Assuming gene expression count matrix format. "
        f"Routing to RNA-seq pipeline: CPM normalisation, "
        f"differential expression analysis, volcano plot, PCA, heatmap."
    )

    decision = RoutingDecision(
        file_type=detection.file_type,
        pipeline_name="RNA-seq Differential Expression Pipeline",
        confidence=detection.confidence,
        reasoning=reasoning
    )

    logger.info(f"Routing to RNA-seq pipeline.")
    result = run_rnaseq_pipeline(
        path,
        control_label=control_label,
        treatment_label=treatment_label,
        output_dir=output_dir,
        use_rag=use_rag
    )
    return result, decision