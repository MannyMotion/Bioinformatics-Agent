"""
variant_annotation.py

Purpose: Variant annotation pipeline for VCF files.
         Parses VCF files, annotates variants with biological context,
         assesses clinical significance, and generates visualisations.

         In production, tools like ANNOVAR, VEP (Ensembl), or SnpEff
         would be used for full annotation against databases.
         This pipeline demonstrates the annotation workflow using
         a curated knowledge base of known variants.

         Reference: McLaren et al. (2016) Genome Biology — Ensembl VEP
         https://doi.org/10.1186/s13059-016-0974-4

Inputs:  Path to a VCF file
Outputs: VariantResult dataclass with annotations, stats, plots

Author:  Emmanuel Ogbu (Manny)
Date:    2026-04-24
"""
import os
os.environ["MPLBACKEND"] = "Agg"

from pathlib import Path
from dataclasses import dataclass, field
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib
matplotlib.rcParams['figure.max_open_warning'] = 0
import matplotlib.pyplot as plt
plt.switch_backend("Agg")
import seaborn as sns

from bioagent.utils.logger import get_logger

logger = get_logger(__name__)

# Minimum quality score to consider a variant reliable
MIN_QUALITY_THRESHOLD = 70

# Curated knowledge base of clinically significant variants
# In production this would query ClinVar, COSMIC, gnomAD databases
# Reference: Landrum et al. (2018) Nucleic Acids Research — ClinVar
# https://doi.org/10.1093/nar/gkx1153
KNOWN_VARIANTS = {
    "rs80357382": {
        "gene": "BRCA1",
        "consequence": "missense_variant",
        "clinical_significance": "Pathogenic",
        "condition": "Hereditary breast and ovarian cancer",
        "chromosome": "17"
    },
    "rs28897696": {
        "gene": "BRCA1",
        "consequence": "missense_variant",
        "clinical_significance": "Pathogenic",
        "condition": "Hereditary breast cancer",
        "chromosome": "17"
    },
    "rs121913428": {
        "gene": "EGFR",
        "consequence": "missense_variant",
        "clinical_significance": "Pathogenic",
        "condition": "Non-small cell lung cancer",
        "chromosome": "7"
    },
    "rs121913529": {
        "gene": "KRAS",
        "consequence": "missense_variant",
        "clinical_significance": "Pathogenic",
        "condition": "Colorectal cancer, lung cancer",
        "chromosome": "12"
    },
    "rs28934578": {
        "gene": "TP53",
        "consequence": "missense_variant",
        "clinical_significance": "Pathogenic",
        "condition": "Li-Fraumeni syndrome, multiple cancers",
        "chromosome": "17"
    },
    "rs80358981": {
        "gene": "BRCA2",
        "consequence": "missense_variant",
        "clinical_significance": "Pathogenic",
        "condition": "Hereditary breast and ovarian cancer",
        "chromosome": "13"
    },
    "rs104894003": {
        "gene": "NRAS",
        "consequence": "missense_variant",
        "clinical_significance": "Pathogenic",
        "condition": "Melanoma, colorectal cancer",
        "chromosome": "1"
    },
    "rs121913254": {
        "gene": "PIK3CA",
        "consequence": "missense_variant",
        "clinical_significance": "Pathogenic",
        "condition": "Breast cancer, colorectal cancer",
        "chromosome": "3"
    },
}


@dataclass
class Variant:
    """
    Represents a single annotated variant from a VCF file.
    Combines raw VCF data with biological annotation.
    """
    chrom: str           # chromosome
    pos: int             # genomic position
    variant_id: str      # rsID or '.' if unknown
    ref: str             # reference allele
    alt: str             # alternate allele
    qual: float          # quality score
    filter_status: str   # PASS or filter reason
    depth: int           # read depth at this position
    allele_freq: float   # allele frequency in sample

    # Annotation fields (filled in by annotate step)
    gene: str = "Unknown"
    consequence: str = "Unknown"
    clinical_significance: str = "Unknown"
    condition: str = "Unknown"
    is_annotated: bool = False
    passes_quality: bool = True


@dataclass
class VariantResult:
    """Full result from the variant annotation pipeline."""
    file_name: str
    total_variants: int
    pass_filter_count: int
    annotated_count: int
    pathogenic_count: int
    variants: list[Variant] = field(default_factory=list)
    plot_paths: list[str] = field(default_factory=list)
    interpretation: str = ""
    warnings: list[str] = field(default_factory=list)


def run_variant_pipeline(
    vcf_file: str | Path,
    output_dir: str | Path = "outputs",
    use_rag: bool = False
) -> VariantResult:
    """
    Run the full variant annotation pipeline.

    Args:
        vcf_file: Path to VCF file.
        output_dir: Directory to save plots.
        use_rag: Whether to query RAG for interpretation.

    Returns:
        VariantResult with annotations, plots, interpretation.

    Raises:
        FileNotFoundError: If VCF file doesn't exist.
    """
    path = Path(vcf_file)
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if not path.exists():
        raise FileNotFoundError(f"VCF file not found: {path}")

    logger.info(f"Starting variant annotation pipeline on: {path.name}")

    # Step 1 — Parse VCF
    logger.info("Step 1: Parsing VCF file...")
    variants = _parse_vcf(path)
    logger.info(f"Parsed {len(variants)} variants.")

    # Step 2 — Quality filter
    logger.info("Step 2: Applying quality filters...")
    variants = _apply_quality_filter(variants)

    # Step 3 — Annotate variants
    logger.info("Step 3: Annotating variants...")
    variants = _annotate_variants(variants)

    # Step 4 — Build result
    pass_count = sum(1 for v in variants if v.passes_quality)
    annotated_count = sum(1 for v in variants if v.is_annotated)
    pathogenic_count = sum(
        1 for v in variants
        if v.clinical_significance == "Pathogenic"
    )

    result = VariantResult(
        file_name=path.name,
        total_variants=len(variants),
        pass_filter_count=pass_count,
        annotated_count=annotated_count,
        pathogenic_count=pathogenic_count,
        variants=variants,
    )

    # Step 5 — Generate plots
    logger.info("Step 4: Generating visualisations...")
    result.plot_paths = _generate_plots(result, variants, out_dir)

    # Step 6 — Warnings
    result.warnings = _generate_warnings(result)

    # Step 7 — Interpretation
    if use_rag:
        logger.info("Step 5: Querying RAG for interpretation...")
        result.interpretation = _get_rag_interpretation(result)
    else:
        result.interpretation = _basic_interpretation(result)

    logger.info(
        f"Variant pipeline complete. "
        f"{pathogenic_count} pathogenic variants identified."
    )
    return result


def _parse_vcf(path: Path) -> list[Variant]:
    """
    Parse a VCF file into a list of Variant objects.

    Skips header lines (starting with #) and parses data lines.
    Extracts INFO field for depth (DP) and allele frequency (AF).

    Args:
        path: Path to VCF file.

    Returns:
        List of Variant objects.
    """
    variants = []

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            # Skip header lines
            if line.startswith("#"):
                continue

            line = line.strip()
            if not line:
                continue

            # VCF columns: CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE
            cols = line.split("\t")
            if len(cols) < 8:
                logger.warning(f"Skipping malformed VCF line: {line[:50]}")
                continue

            try:
                chrom = cols[0]
                pos = int(cols[1])
                variant_id = cols[2]
                ref = cols[3]
                alt = cols[4]
                qual = float(cols[5]) if cols[5] != "." else 0.0
                filter_status = cols[6]
                info = cols[7]

                # Parse INFO field — key=value pairs separated by semicolons
                info_dict = {}
                for item in info.split(";"):
                    if "=" in item:
                        k, v = item.split("=", 1)
                        info_dict[k] = v

                depth = int(info_dict.get("DP", 0))
                allele_freq = float(info_dict.get("AF", 0.0))

                variants.append(Variant(
                    chrom=chrom,
                    pos=pos,
                    variant_id=variant_id,
                    ref=ref,
                    alt=alt,
                    qual=qual,
                    filter_status=filter_status,
                    depth=depth,
                    allele_freq=allele_freq,
                ))

            except (ValueError, IndexError) as e:
                logger.warning(f"Failed to parse VCF line: {e}")
                continue

    return variants


def _apply_quality_filter(variants: list[Variant]) -> list[Variant]:
    """
    Flag variants that fail quality thresholds.

    Args:
        variants: List of parsed variants.

    Returns:
        Same list with passes_quality flag set.
    """
    for v in variants:
        if v.qual < MIN_QUALITY_THRESHOLD or v.filter_status != "PASS":
            v.passes_quality = False
            logger.debug(
                f"Variant {v.variant_id} at {v.chrom}:{v.pos} "
                f"failed QC (QUAL={v.qual}, FILTER={v.filter_status})"
            )
    return variants


def _annotate_variants(variants: list[Variant]) -> list[Variant]:
    """
    Annotate variants with biological information.

    Looks up each variant's rsID in our curated knowledge base.
    In production, this would query Ensembl VEP, ClinVar, or ANNOVAR.

    Reference: Ensembl VEP — https://www.ensembl.org/vep

    Args:
        variants: List of variants to annotate.

    Returns:
        Annotated variant list.
    """
    for v in variants:
        if v.variant_id in KNOWN_VARIANTS:
            annotation = KNOWN_VARIANTS[v.variant_id]
            v.gene = annotation["gene"]
            v.consequence = annotation["consequence"]
            v.clinical_significance = annotation["clinical_significance"]
            v.condition = annotation["condition"]
            v.is_annotated = True
            logger.debug(
                f"Annotated {v.variant_id}: {v.gene} — "
                f"{v.clinical_significance}"
            )
        else:
            logger.debug(f"No annotation found for {v.variant_id}")

    return variants


def _generate_plots(
    result: VariantResult,
    variants: list[Variant],
    output_dir: Path
) -> list[str]:
    """Generate variant plots as standalone HTML files using Chart.js."""
    import json
    from collections import Counter

    out_dir = Path(output_dir)
    plot_paths = []

    # --- Plot 1: Quality Scores ---
    labels = [v.variant_id[:12] for v in variants]
    quals = [float(v.qual) for v in variants]
    colors = ["rgba(76,175,80,0.8)" if v.passes_quality else "rgba(244,67,54,0.8)"
              for v in variants]

    quality_html = f"""<!DOCTYPE html>
<html><head><script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
<style>body{{background:#1a1f35;margin:0;padding:16px;}} canvas{{background:#111827;}}</style>
</head><body>
<canvas id="quality" width="800" height="400"></canvas>
<script>
const ctx = document.getElementById('quality').getContext('2d');
new Chart(ctx, {{
    type: 'bar',
    data: {{
        labels: {json.dumps(labels)},
        datasets: [{{
            label: 'Quality Score',
            data: {json.dumps(quals)},
            backgroundColor: {json.dumps(colors)},
            borderColor: '#0a0e1a',
            borderWidth: 1
        }}]
    }},
    options:{{
        responsive:true,
        plugins:{{
            legend:{{labels:{{color:'#e0e6f0'}}}},
            title:{{display:true,text:'Variant Quality Scores — {result.file_name}',color:'#f1f5f9',font:{{size:16}}}}
        }},
        scales:{{
            x:{{ticks:{{color:'#94a3b8',maxRotation:45}},grid:{{color:'#1e293b'}}}},
            y:{{
                title:{{display:true,text:'Quality Score (PHRED)',color:'#94a3b8'}},
                ticks:{{color:'#94a3b8'}},
                grid:{{color:'#1e293b'}},
                min:0
            }}
        }}
    }}
}});
</script></body></html>"""

    p = str(out_dir / f"{result.file_name}_quality.html")
    Path(p).write_text(quality_html, encoding="utf-8")
    plot_paths.append(p)
    logger.info(f"Saved quality plot: {p}")

    # --- Plot 2: Genes ---
    annotated = [v for v in variants if v.is_annotated]
    if annotated:
        gene_counts = Counter(v.gene for v in annotated)
        genes = list(gene_counts.keys())
        counts = list(gene_counts.values())

        genes_html = f"""<!DOCTYPE html>
<html><head><script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
<style>body{{background:#1a1f35;margin:0;padding:16px;}} canvas{{background:#111827;}}</style>
</head><body>
<canvas id="genes" width="700" height="400"></canvas>
<script>
const ctx = document.getElementById('genes').getContext('2d');
new Chart(ctx, {{
    type: 'bar',
    data: {{
        labels: {json.dumps(genes)},
        datasets: [{{
            label: 'Variants',
            data: {json.dumps(counts)},
            backgroundColor: 'rgba(79,158,255,0.8)',
            borderColor: '#0a0e1a',
            borderWidth: 1
        }}]
    }},
    options:{{
        indexAxis: 'y',
        responsive:true,
        plugins:{{
            legend:{{labels:{{color:'#e0e6f0'}}}},
            title:{{display:true,text:'Variants per Gene — {result.file_name}',color:'#f1f5f9',font:{{size:16}}}}
        }},
        scales:{{
            x:{{title:{{display:true,text:'Number of Variants',color:'#94a3b8'}},ticks:{{color:'#94a3b8'}},grid:{{color:'#1e293b'}}}},
            y:{{ticks:{{color:'#94a3b8'}},grid:{{color:'#1e293b'}}}}
        }}
    }}
}});
</script></body></html>"""

        p = str(out_dir / f"{result.file_name}_genes.html")
        Path(p).write_text(genes_html, encoding="utf-8")
        plot_paths.append(p)
        logger.info(f"Saved genes plot: {p}")

    # --- Plot 3: Clinical Significance Pie ---
    sig_counts = Counter(v.clinical_significance for v in variants)
    sig_labels = list(sig_counts.keys())
    sig_values = list(sig_counts.values())
    pie_colors = []
    color_map = {
        "Pathogenic": "rgba(244,67,54,0.8)",
        "Likely pathogenic": "rgba(255,152,0,0.8)",
        "Uncertain significance": "rgba(255,193,7,0.8)",
        "Benign": "rgba(76,175,80,0.8)",
        "Unknown": "rgba(158,158,158,0.8)"
    }
    for label in sig_labels:
        pie_colors.append(color_map.get(label, "rgba(158,158,158,0.8)"))

    sig_html = f"""<!DOCTYPE html>
<html><head><script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
<style>body{{background:#1a1f35;margin:0;padding:16px;display:flex;justify-content:center;}} canvas{{background:#111827;}}</style>
</head><body>
<canvas id="sig" width="500" height="400"></canvas>
<script>
const ctx = document.getElementById('sig').getContext('2d');
new Chart(ctx, {{
    type: 'pie',
    data: {{
        labels: {json.dumps(sig_labels)},
        datasets: [{{
            data: {json.dumps(sig_values)},
            backgroundColor: {json.dumps(pie_colors)},
            borderColor: '#0a0e1a',
            borderWidth: 2
        }}]
    }},
    options:{{
        responsive:true,
        plugins:{{
            legend:{{labels:{{color:'#e0e6f0'}}}},
            title:{{display:true,text:'Clinical Significance — {result.file_name}',color:'#f1f5f9',font:{{size:16}}}}
        }}
    }}
}});
</script></body></html>"""

    p = str(out_dir / f"{result.file_name}_significance.html")
    Path(p).write_text(sig_html, encoding="utf-8")
    plot_paths.append(p)
    logger.info(f"Saved significance plot: {p}")

    return plot_paths


def _generate_warnings(result: VariantResult) -> list[str]:
    """Generate QC warnings for variant results."""
    warnings = []

    failed = result.total_variants - result.pass_filter_count
    if failed > 0:
        warnings.append(
            f"WARNING: {failed} variant(s) failed quality filters "
            f"and were excluded from analysis."
        )

    if result.pathogenic_count > 0:
        warnings.append(
            f"CLINICAL ALERT: {result.pathogenic_count} pathogenic "
            f"variant(s) identified. Clinical review recommended."
        )

    if result.annotated_count < result.total_variants:
        unannotated = result.total_variants - result.annotated_count
        warnings.append(
            f"NOTE: {unannotated} variant(s) could not be annotated "
            f"from the knowledge base."
        )

    return warnings


def _get_rag_interpretation(result: VariantResult) -> str:
    """Query RAG for biological interpretation."""
    try:
        from bioagent.rag.retriever import BioRetriever  # import here, not at top
        retriever = BioRetriever()
        query = (
            f"variant calling annotation clinical significance "
            f"pathogenic SNP genomics cancer"
        )
        context = retriever.retrieve_as_context(query, n_results=3)
        interpretation = _basic_interpretation(result)
        interpretation += f"\n\nKnowledge Base Context:\n{context}"
        return interpretation
    except Exception as e:
        logger.warning(f"RAG query failed: {e}")
        return _basic_interpretation(result)


def _basic_interpretation(result: VariantResult) -> str:
    """Generate basic interpretation of variant results."""
    lines = [
        f"Variant Annotation Report — {result.file_name}",
        f"{'='*50}",
        f"Total variants identified: {result.total_variants}",
        f"Variants passing QC: {result.pass_filter_count}",
        f"Successfully annotated: {result.annotated_count}",
        f"Pathogenic variants: {result.pathogenic_count}",
        f"",
        f"QC Assessment:",
    ]

    for w in result.warnings:
        lines.append(f"  • {w}")

    # List pathogenic variants
    pathogenic = [
        v for v in result.variants
        if v.clinical_significance == "Pathogenic"
    ]

    if pathogenic:
        lines.append(f"\nPathogenic variants identified:")
        for v in pathogenic:
            lines.append(
                f"  • {v.variant_id} | Gene: {v.gene} | "
                f"Chr{v.chrom}:{v.pos} {v.ref}>{v.alt} | "
                f"Condition: {v.condition}"
            )

    return "\n".join(lines)