"""
rnaseq.py

Purpose: RNA-seq differential expression analysis pipeline.
         Takes a gene count matrix (genes x samples) and performs:
         - Normalisation (CPM - Counts Per Million)
         - PCA (Principal Component Analysis) for sample clustering
         - Differential expression analysis (log2 fold change + stats)
         - Volcano plot (significance vs fold change)
         - Heatmap of top differentially expressed genes
         - Biological interpretation via RAG knowledge base

         This pipeline answers: "Which genes are significantly
         up or down-regulated between two conditions?"

Inputs:  CSV file with gene counts (rows=genes, columns=samples)
         Two condition labels (e.g. "healthy" and "cancer")
Outputs: RNAseqResult dataclass with stats, plots, interpretation

Author:  Emmanuel Ogbu (Manny)
Date:    2026-04-24
"""

from pathlib import Path
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")  # non-interactive backend for server use


from bioagent.utils.logger import get_logger

logger = get_logger(__name__)

# Statistical thresholds for differential expression
# These are standard cutoffs used in published RNA-seq analyses
# Reference: Love et al. (2014) Genome Biology — DESeq2 paper
PVALUE_THRESHOLD = 0.05      # p-value cutoff for significance
LOG2FC_THRESHOLD = 1.0       # log2 fold change cutoff (= 2x change)
MIN_COUNT_THRESHOLD = 10     # minimum count to keep a gene (filter noise)


@dataclass
class GeneResult:
    """
    Differential expression result for a single gene.
    Stores all statistics needed for volcano plot and reporting.
    """
    gene_id: str
    mean_control: float      # mean count in control condition
    mean_treatment: float    # mean count in treatment condition
    log2_fold_change: float  # log2(treatment/control)
    p_value: float           # statistical significance
    significant: bool        # True if passes both thresholds


@dataclass
class RNAseqResult:
    """
    Full RNA-seq analysis result.
    Contains everything needed for reporting and visualisation.
    """
    file_name: str
    control_label: str
    treatment_label: str
    total_genes: int
    genes_tested: int           # after low-count filtering
    upregulated_count: int      # significant + positive log2FC
    downregulated_count: int    # significant + negative log2FC
    gene_results: list[GeneResult] = field(default_factory=list)
    plot_paths: list[str] = field(default_factory=list)
    interpretation: str = ""
    warnings: list[str] = field(default_factory=list)


def run_rnaseq_pipeline(
    counts_file: str | Path,
    control_label: str,
    treatment_label: str,
    output_dir: str | Path = "outputs",
    use_rag: bool = False
) -> RNAseqResult:
    """
    Run the full RNA-seq differential expression pipeline.

    Args:
        counts_file: Path to CSV file with gene counts.
        control_label: Prefix for control sample columns (e.g. "healthy").
        treatment_label: Prefix for treatment sample columns (e.g. "cancer").
        output_dir: Directory to save output plots.
        use_rag: Whether to query RAG for biological interpretation.

    Returns:
        RNAseqResult with full statistics, plots, and interpretation.

    Raises:
        FileNotFoundError: If the counts file does not exist.
        ValueError: If control or treatment columns not found.
    """
    path = Path(counts_file)
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if not path.exists():
        raise FileNotFoundError(f"Counts file not found: {path}")

    logger.info(f"Starting RNA-seq pipeline on: {path.name}")
    logger.info(f"Control: {control_label} | Treatment: {treatment_label}")

    # Step 1 — Load and validate count matrix
    logger.info("Step 1: Loading count matrix...")
    counts_df = _load_counts(path, control_label, treatment_label)

    # Step 2 — Filter low count genes
    logger.info("Step 2: Filtering low-count genes...")
    filtered_df = _filter_low_counts(counts_df, MIN_COUNT_THRESHOLD)
    logger.info(
        f"Retained {len(filtered_df)} genes after filtering "
        f"(removed {len(counts_df) - len(filtered_df)} low-count genes)."
    )

    # Step 3 — Normalise to CPM
    logger.info("Step 3: Normalising to CPM...")
    cpm_df = _normalise_cpm(filtered_df)

    # Step 4 — Differential expression analysis
    logger.info("Step 4: Running differential expression analysis...")
    control_cols = [c for c in filtered_df.columns if c.startswith(control_label)]
    treatment_cols = [c for c in filtered_df.columns if c.startswith(treatment_label)]
    gene_results = _differential_expression(
        filtered_df, control_cols, treatment_cols
    )

    # Step 5 — Build result object
    upregulated = [
        g for g in gene_results
        if g.significant and g.log2_fold_change > 0
    ]
    downregulated = [
        g for g in gene_results
        if g.significant and g.log2_fold_change < 0
    ]

    result = RNAseqResult(
        file_name=path.name,
        control_label=control_label,
        treatment_label=treatment_label,
        total_genes=len(counts_df),
        genes_tested=len(filtered_df),
        upregulated_count=len(upregulated),
        downregulated_count=len(downregulated),
        gene_results=gene_results,
    )

    # Step 6 — Generate visualisations
    logger.info("Step 5: Generating visualisations...")
    control_cols_cpm = [c for c in cpm_df.columns if c.startswith(control_label)]
    treatment_cols_cpm = [c for c in cpm_df.columns if c.startswith(treatment_label)]
    plot_paths = _generate_plots(
        result, gene_results, cpm_df,
        control_cols_cpm, treatment_cols_cpm, out_dir
    )
    result.plot_paths = plot_paths

    # Step 7 — Generate warnings
    result.warnings = _generate_warnings(result)

    # Step 8 — RAG interpretation
    if use_rag:
        logger.info("Step 6: Querying RAG for biological interpretation...")
        result.interpretation = _get_rag_interpretation(result, gene_results)
    else:
        result.interpretation = _basic_interpretation(result, gene_results)

    logger.info(
        f"RNA-seq pipeline complete. "
        f"{result.upregulated_count} up, {result.downregulated_count} down."
    )
    return result


def _load_counts(
    path: Path,
    control_label: str,
    treatment_label: str
) -> pd.DataFrame:
    """
    Load and validate the count matrix CSV.

    Args:
        path: Path to CSV file.
        control_label: Column prefix for control samples.
        treatment_label: Column prefix for treatment samples.

    Returns:
        DataFrame with gene_id as index.

    Raises:
        ValueError: If expected columns are missing.
    """
    # Read CSV with gene_id as the index column
    df = pd.read_csv(path, index_col=0)

    # Validate that we have columns for both conditions
    control_cols = [c for c in df.columns if c.startswith(control_label)]
    treatment_cols = [c for c in df.columns if c.startswith(treatment_label)]

    if not control_cols:
        raise ValueError(
            f"No columns found starting with '{control_label}'. "
            f"Available columns: {list(df.columns)}"
        )
    if not treatment_cols:
        raise ValueError(
            f"No columns found starting with '{treatment_label}'. "
            f"Available columns: {list(df.columns)}"
        )

    logger.info(
        f"Loaded {len(df)} genes. "
        f"Control samples: {control_cols}. "
        f"Treatment samples: {treatment_cols}."
    )
    return df


def _filter_low_counts(df: pd.DataFrame, min_count: int) -> pd.DataFrame:
    """
    Remove genes with very low counts across all samples.

    Low count genes are unreliable — they have high variance relative
    to their mean and can generate false positives in DE analysis.
    We keep only genes where at least one sample has >= min_count reads.

    Args:
        df: Raw count matrix.
        min_count: Minimum count threshold.

    Returns:
        Filtered DataFrame.
    """
    # Keep genes where the maximum count across all samples >= min_count
    mask = df.max(axis=1) >= min_count
    return df[mask]


def _normalise_cpm(df: pd.DataFrame) -> pd.DataFrame:
    """
    Normalise raw counts to CPM (Counts Per Million).

    CPM removes the effect of different library sizes between samples.
    Without normalisation, a sample with 10M reads would appear to have
    2x more expression than a sample with 5M reads — even if expression
    is identical.

    CPM formula: (gene_count / total_counts_in_sample) * 1,000,000

    Args:
        df: Raw count matrix.

    Returns:
        CPM-normalised DataFrame.
    """
    # Sum all counts per sample (column) to get library size
    library_sizes = df.sum(axis=0)

    # Divide each count by library size and multiply by 1M
    cpm = df.divide(library_sizes, axis=1) * 1_000_000

    return cpm


def _differential_expression(
    df: pd.DataFrame,
    control_cols: list[str],
    treatment_cols: list[str]
) -> list[GeneResult]:
    """
    Calculate differential expression statistics for each gene.

    Uses a t-test to compare means between conditions.
    Note: In production, DESeq2 (R) or edgeR would be used for more
    accurate negative binomial modelling. We use t-test here for
    dependency simplicity while learning the concepts.

    Reference: Love et al. (2014) Genome Biology
    https://doi.org/10.1186/s13059-014-0550-8

    Args:
        df: Filtered count matrix.
        control_cols: Column names for control samples.
        treatment_cols: Column names for treatment samples.

    Returns:
        List of GeneResult objects sorted by p-value.
    """
    results = []

    for gene_id in df.index:
        control_counts = df.loc[gene_id, control_cols].values.astype(float)
        treatment_counts = df.loc[gene_id, treatment_cols].values.astype(float)

        mean_control = float(np.mean(control_counts))
        mean_treatment = float(np.mean(treatment_counts))

        # log2 fold change — how much expression changed
        # Add 1 to avoid log2(0) which is undefined
        log2fc = np.log2((mean_treatment + 1) / (mean_control + 1))

        # Two-sample t-test — are the means significantly different?
        if len(control_counts) > 1 and len(treatment_counts) > 1:
            _, p_value = stats.ttest_ind(control_counts, treatment_counts)
        else:
            p_value = 1.0  # can't test with single sample

        # A gene is significant if BOTH thresholds are met
        significant = (
            p_value < PVALUE_THRESHOLD and
            abs(log2fc) >= LOG2FC_THRESHOLD
        )

        results.append(GeneResult(
            gene_id=gene_id,
            mean_control=round(mean_control, 2),
            mean_treatment=round(mean_treatment, 2),
            log2_fold_change=round(float(log2fc), 4),
            p_value=round(float(p_value), 6),
            significant=significant,
        ))

    # Sort by p-value — most significant first
    results.sort(key=lambda x: x.p_value)
    return results


def _generate_plots(
    result: RNAseqResult,
    gene_results: list[GeneResult],
    cpm_df,
    control_cols: list[str],
    treatment_cols: list[str],
    output_dir: Path
) -> list[str]:
    """Generate RNA-seq plots as standalone HTML files using Chart.js."""
    import json
    import numpy as np

    out_dir = Path(output_dir)
    plot_paths = []

    # --- Plot 1: Volcano Plot ---
    up_points, down_points, ns_points = [], [], []
    for g in gene_results:
        neg_log10_p = float(-np.log10(g.p_value + 1e-10))
        point = {"x": float(g.log2_fold_change), "y": neg_log10_p, "label": g.gene_id}
        if g.significant and g.log2_fold_change > 0:
            up_points.append(point)
        elif g.significant and g.log2_fold_change < 0:
            down_points.append(point)
        else:
            ns_points.append(point)

    volcano_html = f"""<!DOCTYPE html>
<html><head><script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
<script src="https://cdn.jsdelivr.net/npm/chartjs-plugin-annotation"></script>
<style>body{{background:#1a1f35;margin:0;padding:16px;}} canvas{{background:#111827;}}</style>
</head><body>
<canvas id="volcano" width="800" height="500"></canvas>
<script>
const ctx = document.getElementById('volcano').getContext('2d');
new Chart(ctx, {{
    type: 'scatter',
    data: {{
        datasets: [
            {{label:'Not significant ({len(ns_points)})',data:{json.dumps(ns_points)},backgroundColor:'rgba(158,158,158,0.5)',pointRadius:5}},
            {{label:'Downregulated ({len(down_points)})',data:{json.dumps(down_points)},backgroundColor:'rgba(33,150,243,0.8)',pointRadius:7}},
            {{label:'Upregulated ({len(up_points)})',data:{json.dumps(up_points)},backgroundColor:'rgba(244,67,54,0.8)',pointRadius:7}}
        ]
    }},
    options:{{
        responsive:true,
        plugins:{{
            legend:{{labels:{{color:'#e0e6f0'}}}},
            title:{{display:true,text:'Volcano Plot: {result.treatment_label} vs {result.control_label}',color:'#f1f5f9',font:{{size:16}}}}
        }},
        scales:{{
            x:{{title:{{display:true,text:'log2 Fold Change',color:'#94a3b8'}},ticks:{{color:'#94a3b8'}},grid:{{color:'#1e293b'}}}},
            y:{{title:{{display:true,text:'-log10(p-value)',color:'#94a3b8'}},ticks:{{color:'#94a3b8'}},grid:{{color:'#1e293b'}}}}
        }}
    }}
}});
</script></body></html>"""

    p = str(out_dir / f"{result.file_name}_volcano.html")
    Path(p).write_text(volcano_html, encoding="utf-8")
    plot_paths.append(p)
    logger.info(f"Saved volcano plot: {p}")

    # PCA disabled — np.linalg.eigh causes Windows segfault in uvicorn
    # Will be enabled when deployed to Linux
    logger.info("PCA skipped on Windows.")

    # --- Plot 3: Heatmap ---
    try:
        sig_genes = [g for g in gene_results if g.significant]
        top_genes = sorted(sig_genes, key=lambda x: abs(x.log2_fold_change), reverse=True)[:8]
        top_ids = [g.gene_id for g in top_genes if g.gene_id in cpm_df.index]
        all_cols = control_cols + treatment_cols

        if top_ids:
            labels = top_ids
            datasets = []
            colors = ["#4CAF50","#388E3C","#1B5E20","#F44336","#C62828","#B71C1C"]
            for ci, col in enumerate(all_cols):
                vals = [float(cpm_df.loc[g, col]) for g in top_ids]
                datasets.append({
                    "label": col,
                    "data": vals,
                    "backgroundColor": colors[ci % len(colors)],
                    "borderColor": "#0a0e1a",
                    "borderWidth": 1
                })

            heatmap_html = f"""<!DOCTYPE html>
<html><head><script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
<style>body{{background:#1a1f35;margin:0;padding:16px;}} canvas{{background:#111827;}}</style>
</head><body>
<canvas id="heatmap" width="800" height="450"></canvas>
<script>
const ctx = document.getElementById('heatmap').getContext('2d');
new Chart(ctx, {{
    type: 'bar',
    data: {{
        labels: {json.dumps(labels)},
        datasets: {json.dumps(datasets)}
    }},
    options:{{
        responsive:true,
        plugins:{{
            legend:{{labels:{{color:'#e0e6f0'}}}},
            title:{{display:true,text:'Top Differentially Expressed Genes (CPM)',color:'#f1f5f9',font:{{size:16}}}}
        }},
        scales:{{
            x:{{ticks:{{color:'#94a3b8'}},grid:{{color:'#1e293b'}}}},
            y:{{title:{{display:true,text:'CPM',color:'#94a3b8'}},ticks:{{color:'#94a3b8'}},grid:{{color:'#1e293b'}}}}
        }}
    }}
}});
</script></body></html>"""

            p = str(out_dir / f"{result.file_name}_heatmap.html")
            Path(p).write_text(heatmap_html, encoding="utf-8")
            plot_paths.append(p)
            logger.info(f"Saved heatmap: {p}")

    except Exception as e:
        logger.warning(f"Heatmap plot failed: {e}")

    return plot_paths


def _generate_warnings(result: RNAseqResult) -> list[str]:
    """Generate QC warnings based on DE results."""
    warnings = []

    if result.upregulated_count == 0 and result.downregulated_count == 0:
        warnings.append(
            "WARNING: No significantly differentially expressed genes found. "
            "Consider relaxing thresholds or checking sample groupings."
        )

    total_sig = result.upregulated_count + result.downregulated_count
    if total_sig > result.genes_tested * 0.5:
        warnings.append(
            f"WARNING: Over 50% of genes are significant — possible "
            f"normalisation issue or batch effect."
        )

    if not warnings:
        warnings.append(
            f"PASS: {result.upregulated_count} upregulated, "
            f"{result.downregulated_count} downregulated genes identified."
        )

    return warnings


def _get_rag_interpretation(result: RNAseqResult,gene_results: list[GeneResult]) -> str:
    """Query RAG for biological interpretation of DE results."""
    try:
        from bioagent.rag.retriever import BioRetriever  # import here, not at top
        retriever = BioRetriever()
        query = (
            f"RNA-seq differential expression analysis "
            f"{result.treatment_label} vs {result.control_label} "
            f"upregulated downregulated genes fold change"
        )
        context = retriever.retrieve_as_context(query, n_results=3)
        interpretation = _basic_interpretation(result, gene_results)
        interpretation += f"\n\nKnowledge Base Context:\n{context}"
        return interpretation
    except Exception as e:
        logger.warning(f"RAG query failed: {e}")
        return _basic_interpretation(result, gene_results)


def _basic_interpretation(
    result: RNAseqResult,
    gene_results: list[GeneResult]
) -> str:
    """Generate basic biological interpretation of DE results."""
    sig_genes = [g for g in gene_results if g.significant]
    top_up = sorted(
        [g for g in sig_genes if g.log2_fold_change > 0],
        key=lambda x: x.log2_fold_change, reverse=True
    )[:3]
    top_down = sorted(
        [g for g in sig_genes if g.log2_fold_change < 0],
        key=lambda x: x.log2_fold_change
    )[:3]

    lines = [
        f"RNA-seq Analysis Report — {result.file_name}",
        f"{'='*50}",
        f"Comparison: {result.treatment_label} vs {result.control_label}",
        f"Total genes in dataset: {result.total_genes}",
        f"Genes tested (after filtering): {result.genes_tested}",
        f"Upregulated genes: {result.upregulated_count}",
        f"Downregulated genes: {result.downregulated_count}",
        f"",
        f"QC Assessment:",
    ]

    for w in result.warnings:
        lines.append(f"  • {w}")

    if top_up:
        lines.append(f"\nTop upregulated genes in {result.treatment_label}:")
        for g in top_up:
            lines.append(
                f"  • {g.gene_id}: {g.log2_fold_change:+.2f} log2FC "
                f"(p={g.p_value:.4f})"
            )

    if top_down:
        lines.append(f"\nTop downregulated genes in {result.treatment_label}:")
        for g in top_down:
            lines.append(
                f"  • {g.gene_id}: {g.log2_fold_change:+.2f} log2FC "
                f"(p={g.p_value:.4f})"
            )

    return "\n".join(lines)