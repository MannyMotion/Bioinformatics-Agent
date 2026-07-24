"""
gene_ids.py

Purpose: Recognise what kind of gene identifiers a dataset uses and, where
         possible, translate Ensembl gene IDs to human-readable symbols.

         Real expression matrices (e.g. from GEO) usually label rows with
         Ensembl IDs like "ENSG00000141510" rather than symbols like "TP53".
         An LLM interpretation is far more useful when it can say "TP53" than
         "ENSG00000141510", so we detect the ID type and map the well-known
         genes we care about.

         This is intentionally a small curated map (no heavy external
         annotation dependency). Unknown IDs pass through unchanged.

Author:  Emmanuel Ogbu (Manny)
Date:    2026-07-02
"""

import re

# Curated Ensembl gene ID -> HGNC symbol map for genes that commonly matter
# in cancer / QC interpretation. Extend as needed.
ENSEMBL_TO_SYMBOL = {
    "ENSG00000141510": "TP53",
    "ENSG00000012048": "BRCA1",
    "ENSG00000139618": "BRCA2",
    "ENSG00000141736": "ERBB2",   # HER2
    "ENSG00000146648": "EGFR",
    "ENSG00000133703": "KRAS",
    "ENSG00000213281": "NRAS",
    "ENSG00000121879": "PIK3CA",
    "ENSG00000091831": "ESR1",    # estrogen receptor
    "ENSG00000148773": "MKI67",   # proliferation marker
    "ENSG00000129514": "FOXA1",
    "ENSG00000171862": "PTEN",
    "ENSG00000171791": "BCL2",
    "ENSG00000136997": "MYC",
    "ENSG00000039068": "CDH1",
    "ENSG00000105329": "TGFB1",
    "ENSG00000232810": "TNF",
    "ENSG00000136244": "IL6",
    "ENSG00000075624": "ACTB",    # housekeeping
    "ENSG00000111640": "GAPDH",   # housekeeping
}

# Regexes for common identifier styles.
_ENSEMBL_RE = re.compile(r"^ENSG\d{11}(\.\d+)?$", re.IGNORECASE)
_SYMBOL_RE = re.compile(r"^[A-Z][A-Z0-9\-]{1,14}$")


def detect_id_type(ids: list[str]) -> str:
    """
    Classify the identifier style used by a list of row IDs.

    Args:
        ids: A sample of row identifiers (e.g. the first ~50 gene IDs).

    Returns:
        One of: "Ensembl gene ID", "gene symbol", "mixed/unknown".
    """
    if not ids:
        return "mixed/unknown"

    sample = [str(i).strip() for i in ids[:50] if str(i).strip()]
    if not sample:
        return "mixed/unknown"

    ensembl = sum(1 for i in sample if _ENSEMBL_RE.match(i))
    symbol = sum(1 for i in sample if _SYMBOL_RE.match(i))

    if ensembl >= 0.6 * len(sample):
        return "Ensembl gene ID"
    if symbol >= 0.6 * len(sample):
        return "gene symbol"
    return "mixed/unknown"


def to_symbol(gene_id: str) -> str:
    """
    Translate a single Ensembl gene ID to its symbol if known.

    Strips any Ensembl version suffix (".3") before lookup. Returns the
    original ID unchanged if it is not a known Ensembl gene.

    Args:
        gene_id: A gene identifier.

    Returns:
        The gene symbol if known, else the original identifier.
    """
    key = str(gene_id).strip().upper().split(".")[0]
    return ENSEMBL_TO_SYMBOL.get(key, gene_id)
