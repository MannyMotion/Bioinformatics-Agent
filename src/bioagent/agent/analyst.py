"""
analyst.py

Purpose: The reasoning layer that turns BioAgent from a fixed router into an
         agent. Given an uploaded file it:

           1. PROFILES the data   — reads it and describes what it is
                                     (dimensions, ID type, candidate sample
                                     groups, cleanliness notes).
           2. PROPOSES questions  — data-driven biological questions the file
                                     can actually answer ("ask a question by
                                     itself"), so the user can pick one or
                                     type their own.
           3. INFERS groups       — works out which columns are control vs
                                     treatment from their names, instead of
                                     forcing the user to use fixed prefixes.
           4. PLANS the analysis  — maps the chosen question + profile to a
                                     concrete pipeline and parameters.

         Profiling and group inference are deterministic and instant, so the
         upload step stays fast and works offline. The heavier LLM work
         (interpretation) still happens later, during analysis.

Author:  Emmanuel Ogbu (Manny)
Date:    2026-07-02
"""

import re
from pathlib import Path
from dataclasses import dataclass, field
from typing import Any

import pandas as pd

from bioagent.agent.detector import DetectionResult
from bioagent.agent import gene_ids
from bioagent.utils.logger import get_logger

logger = get_logger(__name__)


# --- Keyword vocabulary for sample-group inference ---
# Matched against whole tokens (split on non-alphanumeric), so "wt" matches
# the "WT" in "8-Breast-WT" but not the "wt" inside "growth".
CONTROL_KEYWORDS = {
    "control", "ctrl", "untreated", "wt", "wildtype", "normal", "healthy",
    "hc", "mock", "ref", "reference", "baseline", "vehicle", "naive", "pre",
}
TREATMENT_KEYWORDS = {
    "treatment", "treated", "treat", "case", "mutant", "mut", "tumor",
    "tumour", "cancer", "disease", "diseased", "ko", "knockout", "patient",
    "stimulated", "infected", "exposed", "post", "drug",
}


@dataclass
class DataProfile:
    """
    A structured description of an uploaded dataset, produced before any
    heavy analysis runs. This is what lets the agent reason about the data.
    """
    file_type: str
    summary: str                       # one-paragraph plain-English overview
    suggested_questions: list[str] = field(default_factory=list)
    groups: dict[str, list[str]] = field(default_factory=dict)  # name -> columns
    group_method: str = ""             # how groups were inferred
    needs_group_confirmation: bool = False
    id_type: str = ""                  # e.g. "Ensembl gene ID"
    details: dict = field(default_factory=dict)  # type-specific extras

    def as_dict(self) -> dict:
        """JSON-serialisable form for the API and job store."""
        return {
            "file_type": self.file_type,
            "summary": self.summary,
            "suggested_questions": self.suggested_questions,
            "groups": self.groups,
            "group_method": self.group_method,
            "needs_group_confirmation": self.needs_group_confirmation,
            "id_type": self.id_type,
            "details": self.details,
        }

    @classmethod
    def from_dict(cls, data: dict) -> "DataProfile":
        """Rebuild a DataProfile from its stored dict (across API phases)."""
        return cls(
            file_type=data.get("file_type", "UNKNOWN"),
            summary=data.get("summary", ""),
            suggested_questions=data.get("suggested_questions", []),
            groups=data.get("groups", {}),
            group_method=data.get("group_method", ""),
            needs_group_confirmation=data.get("needs_group_confirmation", False),
            id_type=data.get("id_type", ""),
            details=data.get("details", {}),
        )


# --- Public entry point ---

def profile_file(file_path: str | Path, detection: DetectionResult) -> DataProfile:
    """
    Profile an uploaded file according to its detected type.

    Args:
        file_path: Path to the uploaded file.
        detection: The DetectionResult from the detector.

    Returns:
        A DataProfile describing the data and proposing questions.
    """
    path = Path(file_path)
    ftype = detection.file_type.upper()
    logger.info(f"Profiling {path.name} as {ftype}...")

    if ftype in ("CSV", "TSV"):
        return _profile_expression_matrix(path, detection)
    if ftype in ("FASTA", "FASTQ"):
        return _profile_fasta(path, detection)
    if ftype == "VCF":
        return _profile_vcf(path, detection)

    # Unknown — nothing to profile, but stay graceful.
    return DataProfile(
        file_type=ftype,
        summary="The file type could not be determined, so no analysis "
                "plan could be prepared.",
        needs_group_confirmation=False,
    )


# --- Expression matrix (CSV/TSV) ---

def _profile_expression_matrix(path: Path, detection: DetectionResult) -> DataProfile:
    """Profile a gene x sample count matrix and infer condition groups."""
    delimiter = detection.metadata.get("delimiter")
    # sep=None sniffs if the detector did not record one.
    df = pd.read_csv(
        path, index_col=0,
        sep=delimiter if delimiter else None,
        engine="python", nrows=2000  # a sample is enough to profile
    )

    columns = [str(c) for c in df.columns]
    id_type = gene_ids.detect_id_type([str(i) for i in df.index[:50]])
    groups, method, confident = infer_groups(columns)

    n_genes_shown = len(df)
    summary = (
        f"This looks like a gene-expression count matrix: "
        f"{len(columns)} sample columns and at least {n_genes_shown} genes "
        f"(row IDs are {id_type}). "
    )
    if confident:
        gnames = list(groups.keys())
        summary += (
            f"I inferred two groups from the column names — "
            f"'{gnames[0]}' ({len(groups[gnames[0]])} samples) vs "
            f"'{gnames[1]}' ({len(groups[gnames[1]])} samples) — "
            f"so I can run a differential expression comparison."
        )
    else:
        summary += (
            "I could not confidently split the samples into two groups from "
            "their names, so please confirm which columns are the control "
            "and which are the treatment before I run the comparison."
        )

    profile = DataProfile(
        file_type=detection.file_type,
        summary=summary,
        groups=groups,
        group_method=method,
        needs_group_confirmation=not confident,
        id_type=id_type,
        details={
            "n_columns": len(columns),
            "n_genes_sampled": n_genes_shown,
            "columns": columns,
        },
    )
    profile.suggested_questions = _questions_for_expression(groups, confident)
    return profile


def infer_groups(columns: list[str]) -> tuple[dict[str, list[str]], str, bool]:
    """
    Infer two sample groups (control vs treatment) from column names.

    Strategy, in order:
      1. Keyword match — assign each column by control/treatment vocabulary.
      2. Two-token design — if names share exactly two recurring condition
         tokens (e.g. "...WT..." vs "...Her2-ampl..."), split on those.

    Args:
        columns: Sample column names.

    Returns:
        (groups, method, confident) where groups maps a group label to its
        columns, method describes how it was decided, and confident is True
        only when both groups are non-empty and cover most columns.
    """
    # --- Strategy 1: control/treatment keyword vocabulary ---
    control_cols, treatment_cols = [], []
    for col in columns:
        tokens = set(re.split(r"[^a-z0-9]+", col.lower()))
        if tokens & CONTROL_KEYWORDS:
            control_cols.append(col)
        elif tokens & TREATMENT_KEYWORDS:
            treatment_cols.append(col)

    if control_cols and treatment_cols:
        assigned = len(control_cols) + len(treatment_cols)
        confident = assigned >= max(2, int(0.6 * len(columns)))
        return (
            {"control": control_cols, "treatment": treatment_cols},
            "matched control/treatment keywords in the column names",
            confident,
        )

    # --- Strategy 2: dominant two-token design ---
    token_freq: dict[str, int] = {}
    for col in columns:
        for tok in set(re.split(r"[^a-z0-9]+", col.lower())):
            if tok and not tok.isdigit() and len(tok) > 1:
                token_freq[tok] = token_freq.get(tok, 0) + 1

    # Tokens that appear in several (but not all) columns are candidate labels.
    candidates = sorted(
        (t for t, n in token_freq.items() if 1 < n < len(columns)),
        key=lambda t: token_freq[t], reverse=True,
    )
    if len(candidates) >= 2:
        a, b = candidates[0], candidates[1]
        group_a = [c for c in columns if a in re.split(r"[^a-z0-9]+", c.lower())]
        group_b = [c for c in columns
                   if b in re.split(r"[^a-z0-9]+", c.lower()) and c not in group_a]
        if group_a and group_b:
            return (
                {a: group_a, b: group_b},
                f"split on the two most common condition labels "
                f"('{a}' vs '{b}')",
                False,  # plausible but worth a human glance
            )

    # --- Nothing worked ---
    return ({}, "could not infer groups from column names", False)


def _questions_for_expression(groups: dict, confident: bool) -> list[str]:
    """Data-driven questions for an expression matrix."""
    if confident and len(groups) == 2:
        g = list(groups.keys())
        return [
            f"Which genes are significantly differentially expressed "
            f"between {g[1]} and {g[0]}?",
            f"What are the most strongly upregulated genes in {g[1]}?",
            f"Are there known cancer-related genes among the differentially "
            f"expressed genes?",
        ]
    return [
        "Which genes differ most between my two sample groups?",
        "What are the top upregulated and downregulated genes?",
        "Is there a clear expression difference between conditions?",
    ]


# --- FASTA / FASTQ ---

def _profile_fasta(path: Path, detection: DetectionResult) -> DataProfile:
    """Profile a FASTA file: sequence count, length range, sequence type."""
    n = 0
    total_len = 0
    min_len = None
    max_len = 0
    cur = 0
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(">"):
                if n > 0:
                    total_len += cur
                    min_len = cur if min_len is None else min(min_len, cur)
                    max_len = max(max_len, cur)
                n += 1
                cur = 0
            elif not line.startswith("@") and not line.startswith("+"):
                cur += len(line.strip())
        if cur:  # final record
            total_len += cur
            min_len = cur if min_len is None else min(min_len, cur)
            max_len = max(max_len, cur)

    seq_type = detection.metadata.get("sequence_type", "sequence")
    summary = (
        f"This is a {detection.file_type} file with {n} sequence(s) "
        f"({seq_type}), lengths from {min_len or 0} to {max_len} bp. "
        f"I can run quality control: GC content, length distribution, "
        f"nucleotide composition and low-complexity detection."
    )
    profile = DataProfile(
        file_type=detection.file_type,
        summary=summary,
        details={"n_sequences": n, "min_length": min_len or 0, "max_length": max_len},
        suggested_questions=[
            "Is the sequence quality good enough for downstream analysis?",
            "Is the GC content within the expected biological range?",
            "Are there any low-complexity or problematic sequences?",
        ],
    )
    return profile


# --- VCF ---

def _profile_vcf(path: Path, detection: DetectionResult) -> DataProfile:
    """Profile a VCF file: variant count and presence of rsIDs."""
    n_variants = 0
    n_rsids = 0
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue
            if not line.strip():
                continue
            n_variants += 1
            cols = line.split("\t")
            if len(cols) > 2 and cols[2].startswith("rs"):
                n_rsids += 1
            if n_variants >= 100000:  # safety cap for profiling
                break

    summary = (
        f"This is a VCF file with {n_variants} variant(s), "
        f"{n_rsids} of which carry an rsID. I can filter by quality, "
        f"annotate against the clinical knowledge base, and flag any "
        f"pathogenic findings."
    )
    profile = DataProfile(
        file_type="VCF",
        summary=summary,
        details={"n_variants": n_variants, "n_rsids": n_rsids},
        suggested_questions=[
            "Are there any pathogenic variants in this sample?",
            "Which genes carry variants, and are any clinically significant?",
            "How many variants pass quality filtering?",
        ],
    )
    return profile


# --- Planning ---

@dataclass
class AnalysisPlan:
    """The concrete steps the agent will run to answer a question."""
    pipeline: str                       # "fasta_qc" | "rnaseq" | "variant"
    question: str
    steps: list[str] = field(default_factory=list)
    params: dict = field(default_factory=dict)


def plan_analysis(
    question: str,
    profile: DataProfile,
    groups: dict[str, list[str]] | None = None,
) -> AnalysisPlan:
    """
    Map a question + data profile to a concrete analysis plan.

    The data type determines the pipeline; the question shapes how the
    result is interpreted. For expression data the (possibly user-confirmed)
    group mapping is attached as parameters.

    Args:
        question: The chosen or user-typed question.
        profile: The DataProfile from profile_file().
        groups: Optional group mapping that overrides the inferred one
            (e.g. after the user edits it in the UI).

    Returns:
        An AnalysisPlan.
    """
    ftype = profile.file_type.upper()

    if ftype in ("CSV", "TSV"):
        groups = groups or profile.groups
        labels = list(groups.keys())
        control = labels[0] if labels else "control"
        treatment = labels[1] if len(labels) > 1 else "treatment"
        return AnalysisPlan(
            pipeline="rnaseq",
            question=question,
            steps=[
                "Load and validate the count matrix",
                f"Assign samples to groups ({control} vs {treatment})",
                "Filter low-count genes",
                "Normalise to CPM",
                "Test for differential expression (per gene)",
                "Build volcano plot and heatmap",
                "Interpret the findings for the question asked",
            ],
            params={
                "control_label": control,
                "treatment_label": treatment,
                "control_cols": groups.get(control, []),
                "treatment_cols": groups.get(treatment, []),
            },
        )

    if ftype in ("FASTA", "FASTQ"):
        return AnalysisPlan(
            pipeline="fasta_qc",
            question=question,
            steps=[
                "Parse the sequences",
                "Compute per-sequence GC content and length",
                "Detect low-complexity sequences",
                "Summarise nucleotide composition",
                "Interpret the QC result for the question asked",
            ],
        )

    if ftype == "VCF":
        return AnalysisPlan(
            pipeline="variant",
            question=question,
            steps=[
                "Parse the variants",
                "Apply quality filtering",
                "Annotate against the clinical knowledge base",
                "Flag pathogenic findings",
                "Interpret the result for the question asked",
            ],
        )

    return AnalysisPlan(pipeline="none", question=question,
                        steps=["No pipeline available for this file type."])


def autopick_question(profile: DataProfile) -> str:
    """Return the agent's own best question when the user lets it decide."""
    if profile.suggested_questions:
        return profile.suggested_questions[0]
    return "What are the key findings in this dataset?"
