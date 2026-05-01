"""
fasta_qc.py

Purpose: FASTA file Quality Control pipeline.
         Analyses any FASTA file and produces:
         - Per-sequence GC content
         - Sequence length distribution
         - Low complexity sequence detection
         - Summary statistics
         - Publication-quality visualisations
         - Human-readable biological interpretation

         This pipeline is the first analysis module in the BioAgent system.
         It queries the RAG knowledge base to provide context-aware
         explanations of every QC decision made.

Inputs:  Path to a FASTA file (.fasta, .fa, .fna, .faa)
Outputs: QCResult dataclass containing stats, plots, and interpretation

Author:  Emmanuel Ogbu (Manny)
Date:    2026-04-23
"""
import os
os.environ["MPLBACKEND"] = "Agg"

import re
from pathlib import Path
from dataclasses import dataclass, field

import numpy as np
import matplotlib
matplotlib.use("Agg")  # non-interactive backend — no display needed on server
import matplotlib
matplotlib.rcParams['figure.max_open_warning'] = 0
import matplotlib.pyplot as plt
plt.switch_backend("Agg")
import matplotlib.patches as mpatches
import seaborn as sns

from bioagent.parsers.fasta_parser import parse_fasta, FastaRecord
from bioagent.utils.logger import get_logger

logger = get_logger(__name__)

# GC content thresholds — based on standard bioinformatics QC practice
# Reference: FastQC documentation (Andrews, Babraham Bioinformatics)
GC_LOW_THRESHOLD = 35.0   # below this is unusually AT-rich
GC_HIGH_THRESHOLD = 65.0  # above this is unusually GC-rich

# Low complexity threshold — sequences with <30% unique kmers flagged
LOW_COMPLEXITY_THRESHOLD = 0.30


@dataclass
class SequenceStats:
    """
    QC statistics for a single FASTA sequence.
    Stored per-sequence so we can report individual outliers.
    """
    seq_id: str
    length: int
    gc_content: float        # percentage 0-100
    is_low_complexity: bool  # True if sequence is repetitive/low complexity
    nucleotide_counts: dict  # {"A": n, "C": n, "G": n, "T": n, "N": n}


@dataclass
class QCResult:
    """
    Full QC result for an entire FASTA file.
    This is what the pipeline returns — everything downstream needs.
    """
    file_name: str
    total_sequences: int
    total_bases: int
    mean_gc: float
    median_length: float
    min_length: int
    max_length: int
    low_complexity_count: int
    sequences: list[SequenceStats] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    interpretation: str = ""      # RAG-backed biological explanation
    plot_paths: list[str] = field(default_factory=list)  # saved plot file paths


def run_fasta_qc(
    file_path: str | Path,
    output_dir: str | Path = "outputs",
    use_rag: bool = False
) -> QCResult:
    """
    Run the full FASTA QC pipeline on a given file.

    Args:
        file_path: Path to the FASTA file to analyse.
        output_dir: Directory to save output plots.
        use_rag: Whether to query RAG for biological interpretation.

    Returns:
        QCResult with full statistics, warnings, plots, and interpretation.

    Raises:
        FileNotFoundError: If the FASTA file does not exist.
    """
    path = Path(file_path)
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Starting FASTA QC pipeline on: {path.name}")

    # Step 1 — Parse the FASTA file
    logger.info("Step 1: Parsing FASTA file...")
    records = parse_fasta(path)
    logger.info(f"Parsed {len(records)} sequences.")

    # Step 2 — Calculate per-sequence statistics
    logger.info("Step 2: Calculating per-sequence statistics...")
    sequence_stats = [_analyse_sequence(record) for record in records]

    # Step 3 — Calculate summary statistics
    logger.info("Step 3: Calculating summary statistics...")
    result = _calculate_summary(path.name, sequence_stats)

    # Step 4 — Generate warnings
    logger.info("Step 4: Generating QC warnings...")
    result.warnings = _generate_warnings(result, sequence_stats)

    # Step 5 — Generate visualisations in separate process
    logger.info("Step 5: Generating visualisations...")
    try:
        import subprocess
        import sys
        import json

        # Write result data to temp file for subprocess
        temp_data = {
            "file_name": result.file_name,
            "mean_gc": result.mean_gc,
            "median_length": result.median_length,
            "gc_values": [s.gc_content for s in sequence_stats],
            "lengths": [s.length for s in sequence_stats],
            "nucleotide_counts": {
                k: sum(s.nucleotide_counts[k] for s in sequence_stats)
                for k in ["A", "C", "G", "T", "N"]
            },
            "output_dir": str(out_dir)
        }

        temp_file = out_dir / f"{result.file_name}_plotdata.json"
        temp_file.write_text(json.dumps(temp_data))

        # Run plot generation as separate subprocess
        plot_script = Path(__file__).parent.parent / "utils" / "plot_runner.py"
        logger.info(f"Plot script path: {plot_script} — exists: {plot_script.exists()}")
        proc = subprocess.run(
            [sys.executable, str(plot_script), str(temp_file)],
            capture_output=True, text=True, timeout=60
        )
        logger.warning(f"Plot subprocess failed. returncode={proc.returncode} stderr='{proc.stderr[:300]}' stdout='{proc.stdout[:300]}'")
        temp_file.unlink(missing_ok=True)

        if proc.returncode == 0:
            plot_paths_raw = json.loads(proc.stdout.strip())
            result.plot_paths = plot_paths_raw
            logger.info(f"Step 5: Generated {len(plot_paths_raw)} plots.")
        else:
            logger.warning(f"Plot subprocess failed: {proc.stderr[:200]}")
            result.plot_paths = []

    except Exception as e:
        logger.warning(f"Plot generation failed: {type(e).__name__}: {e}")
        result.plot_paths = []

    # Step 6 — Query RAG for biological interpretation
    if use_rag:
        logger.info("Step 6: Querying RAG knowledge base for interpretation...")
        result.interpretation = _get_rag_interpretation(result)
    else:
        result.interpretation = _basic_interpretation(result)

    logger.info(f"FASTA QC complete. {len(result.warnings)} warning(s) generated.")
    return result


def _analyse_sequence(record: FastaRecord) -> SequenceStats:
    """
    Calculate QC metrics for a single sequence.

    Args:
        record: A FastaRecord object from the parser.

    Returns:
        SequenceStats with GC content, complexity, nucleotide counts.
    """
    seq = record.sequence.upper()

    # Count each nucleotide — N represents unknown bases
    counts = {
        "A": seq.count("A"),
        "C": seq.count("C"),
        "G": seq.count("G"),
        "T": seq.count("T"),
        "N": seq.count("N"),
    }

    # GC content = (G + C) / total bases * 100
    # We exclude N from total to avoid penalising sequences with unknown bases
    known_bases = counts["A"] + counts["C"] + counts["G"] + counts["T"]
    gc_content = 0.0
    if known_bases > 0:
        gc_content = (counts["G"] + counts["C"]) / known_bases * 100

    # Low complexity detection using linguistic complexity
    # A sequence is low complexity if it uses very few unique k-mers
    is_low_complexity = _check_low_complexity(seq)

    return SequenceStats(
        seq_id=record.seq_id,
        length=record.length,
        gc_content=round(gc_content, 2),
        is_low_complexity=is_low_complexity,
        nucleotide_counts=counts,
    )


def _check_low_complexity(sequence: str, k: int = 4) -> bool:
    """
    Detect low complexity sequences using k-mer diversity.

    A low complexity sequence (e.g. ATATATATAT or AAAAAAAAAA) has very
    few unique k-mers relative to its length. This can cause issues in
    alignment and assembly steps downstream.

    Args:
        sequence: The raw sequence string.
        k: K-mer size (default 4).

    Returns:
        True if sequence is low complexity, False otherwise.
    """
    if len(sequence) < k:
        return False

    # Generate all k-mers from the sequence
    kmers = [sequence[i:i+k] for i in range(len(sequence) - k + 1)]
    total_kmers = len(kmers)
    unique_kmers = len(set(kmers))

    # Ratio of unique to total k-mers — low ratio = low complexity
    complexity_ratio = unique_kmers / total_kmers
    return complexity_ratio < LOW_COMPLEXITY_THRESHOLD


def _calculate_summary(
    file_name: str,
    stats: list[SequenceStats]
) -> QCResult:
    """
    Calculate summary statistics across all sequences.

    Args:
        file_name: Name of the source file.
        stats: List of per-sequence SequenceStats.

    Returns:
        QCResult with summary statistics populated.
    """
    lengths = [s.length for s in stats]
    gc_values = [s.gc_content for s in stats]
    low_complexity_count = sum(1 for s in stats if s.is_low_complexity)

    return QCResult(
        file_name=file_name,
        total_sequences=len(stats),
        total_bases=sum(lengths),
        mean_gc=round(float(np.mean(gc_values)), 2),
        median_length=float(np.median(lengths)),
        min_length=min(lengths),
        max_length=max(lengths),
        low_complexity_count=low_complexity_count,
        sequences=stats,
    )


def _generate_warnings(
    result: QCResult,
    stats: list[SequenceStats]
) -> list[str]:
    """
    Generate QC warnings based on analysis results.

    Args:
        result: The QCResult summary.
        stats: Per-sequence statistics.

    Returns:
        List of warning strings.
    """
    warnings = []

    # Warning 1 — Overall GC content outside normal range
    if result.mean_gc < GC_LOW_THRESHOLD:
        warnings.append(
            f"WARNING: Mean GC content ({result.mean_gc}%) is below "
            f"the expected threshold ({GC_LOW_THRESHOLD}%). "
            f"This may indicate AT-rich organisms or contamination."
        )
    elif result.mean_gc > GC_HIGH_THRESHOLD:
        warnings.append(
            f"WARNING: Mean GC content ({result.mean_gc}%) is above "
            f"the expected threshold ({GC_HIGH_THRESHOLD}%). "
            f"This may indicate GC-rich organisms or library bias."
        )

    # Warning 2 — Low complexity sequences detected
    if result.low_complexity_count > 0:
        pct = round(result.low_complexity_count / result.total_sequences * 100, 1)
        warnings.append(
            f"WARNING: {result.low_complexity_count} sequence(s) ({pct}%) "
            f"flagged as low complexity. These may cause issues in "
            f"downstream alignment or assembly."
        )

    # Warning 3 — High sequence length variation
    if result.max_length > 0 and result.min_length > 0:
        length_ratio = result.max_length / result.min_length
        if length_ratio > 10:
            warnings.append(
                f"WARNING: Large variation in sequence lengths detected "
                f"(min: {result.min_length}bp, max: {result.max_length}bp). "
                f"Verify this is expected for your data type."
            )

    # Warning 4 — Individual sequences with extreme GC content
    outliers = [
        s for s in stats
        if s.gc_content < GC_LOW_THRESHOLD or s.gc_content > GC_HIGH_THRESHOLD
    ]
    if outliers:
        warnings.append(
            f"WARNING: {len(outliers)} sequence(s) have extreme GC content. "
            f"First outlier: {outliers[0].seq_id} ({outliers[0].gc_content}% GC)."
        )

    if not warnings:
        warnings.append("PASS: All QC metrics within acceptable ranges.")

    return warnings


def _generate_plots(
    result: QCResult,
    stats: list[SequenceStats],
    output_dir: Path
) -> list[str]:
    """
    Generate publication-quality QC visualisations.

    Produces:
    1. GC content distribution histogram
    2. Sequence length distribution
    3. Nucleotide composition bar chart

    Args:
        result: QCResult summary statistics.
        stats: Per-sequence statistics.
        output_dir: Directory to save plots.

    Returns:
        List of file paths to saved plots.
    """
    # Use seaborn style for publication-quality aesthetics
    sns.set_theme(style="whitegrid", palette="muted")
    plot_paths = []

    # --- Plot 1: GC Content Distribution ---
    fig, ax = plt.subplots(figsize=(10, 6))
    gc_values = [s.gc_content for s in stats]

    ax.hist(gc_values, bins=20, color="#2196F3", edgecolor="white", alpha=0.8)
    ax.axvline(
        result.mean_gc, color="#F44336", linestyle="--",
        linewidth=2, label=f"Mean GC: {result.mean_gc}%"
    )
    ax.axvline(
        GC_LOW_THRESHOLD, color="#FF9800", linestyle=":",
        linewidth=1.5, label=f"Low threshold ({GC_LOW_THRESHOLD}%)"
    )
    ax.axvline(
        GC_HIGH_THRESHOLD, color="#FF9800", linestyle=":",
        linewidth=1.5, label=f"High threshold ({GC_HIGH_THRESHOLD}%)"
    )

    ax.set_xlabel("GC Content (%)", fontsize=12)
    ax.set_ylabel("Number of Sequences", fontsize=12)
    ax.set_title(
        f"GC Content Distribution\n{result.file_name} "
        f"({result.total_sequences} sequences)",
        fontsize=14, fontweight="bold"
    )
    ax.legend(fontsize=10)
    plt.tight_layout()

    gc_plot_path = output_dir / f"{result.file_name}_gc_content.png"
    plt.savefig(gc_plot_path, dpi=150, bbox_inches="tight")
    plt.close()
    plot_paths.append(str(gc_plot_path))
    logger.info(f"Saved GC content plot: {gc_plot_path}")

    # --- Plot 2: Sequence Length Distribution ---
    fig, ax = plt.subplots(figsize=(10, 6))
    lengths = [s.length for s in stats]

    ax.hist(lengths, bins=20, color="#4CAF50", edgecolor="white", alpha=0.8)
    ax.axvline(
        result.median_length, color="#F44336", linestyle="--",
        linewidth=2, label=f"Median: {result.median_length:.0f}bp"
    )

    ax.set_xlabel("Sequence Length (bp)", fontsize=12)
    ax.set_ylabel("Number of Sequences", fontsize=12)
    ax.set_title(
        f"Sequence Length Distribution\n{result.file_name}",
        fontsize=14, fontweight="bold"
    )
    ax.legend(fontsize=10)
    plt.tight_layout()

    length_plot_path = output_dir / f"{result.file_name}_lengths.png"
    plt.savefig(length_plot_path, dpi=150, bbox_inches="tight")
    plt.close()
    plot_paths.append(str(length_plot_path))
    logger.info(f"Saved length distribution plot: {length_plot_path}")

    # --- Plot 3: Nucleotide Composition ---
    fig, ax = plt.subplots(figsize=(10, 6))

    # Sum nucleotide counts across all sequences
    total_counts = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0}
    for s in stats:
        for nuc, count in s.nucleotide_counts.items():
            total_counts[nuc] += count

    total = sum(total_counts.values())
    percentages = {k: v/total*100 for k, v in total_counts.items() if total > 0}

    colors = {
        "A": "#4CAF50", "T": "#2196F3",
        "G": "#FF9800", "C": "#F44336", "N": "#9E9E9E"
    }
    bars = ax.bar(
        percentages.keys(),
        percentages.values(),
        color=[colors[k] for k in percentages.keys()],
        edgecolor="white",
        alpha=0.85
    )

    # Add percentage labels on top of each bar
    for bar, (nuc, pct) in zip(bars, percentages.items()):
        ax.text(
            bar.get_x() + bar.get_width()/2,
            bar.get_height() + 0.5,
            f"{pct:.1f}%",
            ha="center", va="bottom", fontsize=11, fontweight="bold"
        )

    ax.set_xlabel("Nucleotide", fontsize=12)
    ax.set_ylabel("Percentage (%)", fontsize=12)
    ax.set_title(
        f"Nucleotide Composition\n{result.file_name}",
        fontsize=14, fontweight="bold"
    )
    plt.tight_layout()

    comp_plot_path = output_dir / f"{result.file_name}_composition.png"
    plt.savefig(comp_plot_path, dpi=150, bbox_inches="tight")
    plt.close()
    plot_paths.append(str(comp_plot_path))
    logger.info(f"Saved nucleotide composition plot: {comp_plot_path}")

    return plot_paths


def _get_rag_interpretation(result: QCResult) -> str:
    """
    Query the RAG knowledge base to generate a biological interpretation
    of the QC results.

    Args:
        result: The QCResult to interpret.

    Returns:
        Human-readable biological interpretation string.
    """
    try:
        from bioagent.rag.retriever import BioRetriever  # import here, not at top
        retriever = BioRetriever()

        # Build a query based on the actual results
        query = (
            f"FASTA sequence quality control GC content {result.mean_gc:.1f}% "
            f"sequence length {result.median_length:.0f}bp bioinformatics"
        )

        context = retriever.retrieve_as_context(query, n_results=3)

        # Build interpretation from results
        interpretation = _basic_interpretation(result)
        interpretation += f"\n\nKnowledge Base Context:\n{context}"

        return interpretation

    except Exception as e:
        logger.warning(f"RAG query failed, using basic interpretation: {e}")
        return _basic_interpretation(result)


def _basic_interpretation(result: QCResult) -> str:
    """
    Generate a basic biological interpretation without RAG.

    Args:
        result: The QCResult to interpret.

    Returns:
        Human-readable interpretation string.
    """
    lines = [
        f"FASTA QC Report — {result.file_name}",
        f"{'='*50}",
        f"Total sequences analysed: {result.total_sequences}",
        f"Total bases: {result.total_bases:,}",
        f"Mean GC content: {result.mean_gc}%",
        f"Median sequence length: {result.median_length:.0f}bp",
        f"Length range: {result.min_length}bp — {result.max_length}bp",
        f"Low complexity sequences: {result.low_complexity_count}",
        f"",
        f"QC Assessment:",
    ]

    for warning in result.warnings:
        lines.append(f"  • {warning}")

    # GC content biological context
    if result.mean_gc < GC_LOW_THRESHOLD:
        lines.append(
            f"\nBiological Note: GC content of {result.mean_gc}% suggests "
            f"an AT-rich organism (e.g. Plasmodium falciparum ~19% GC) or "
            f"possible contamination. Verify sample origin."
        )
    elif result.mean_gc > GC_HIGH_THRESHOLD:
        lines.append(
            f"\nBiological Note: GC content of {result.mean_gc}% suggests "
            f"a GC-rich organism (e.g. Streptomyces ~72% GC) or "
            f"possible PCR bias. Check library preparation."
        )
    else:
        lines.append(
            f"\nBiological Note: GC content of {result.mean_gc}% is within "
            f"the typical range for most organisms (35-65%). "
            f"Human genome GC content is approximately 41%."
        )

    return "\n".join(lines)