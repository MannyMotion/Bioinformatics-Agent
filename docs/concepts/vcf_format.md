# VCF format (Variant Call Format)

**What it is:** The standard text file listing where a sample's DNA differs from a reference genome.

**Structure:**
- `##` lines: metadata (file format version, reference, definitions of INFO tags).
- One `#CHROM …` header line.
- One line per variant, tab-separated. Fixed columns: **CHROM** (chromosome), **POS** (1-based position), **ID** (rsID or `.`), **REF** (reference base(s)), **ALT** (alternate base(s)), **QUAL** (Phred-scaled confidence a variant exists), **FILTER** (`PASS` or the failed filter's name), **INFO** (`key=value;key=value`, e.g. `DP=100;AF=0.5`), then optional FORMAT and per-sample columns.

**Biology / intuition:** A SNV is one base swapped (A→G). An indel is an insertion or deletion. QUAL 30 means about a 1-in-1,000 chance the call is wrong. DP is read depth, AF is the allele fraction.

**When it's wrong to use / pitfalls:**
- Multiple ALT alleles are comma-separated in one line. Our parser stores the string as-is.
- `INFO` keys vary by caller, so never assume DP/AF exist (ours defaults to 0).
- Our quality filter (QUAL ≥ 70 and FILTER = PASS) is a project choice, not a universal standard.
- Position alone isn't enough to name a variant. You need chromosome, position, ref, alt and the genome build.

**Where in our code:** `_parse_vcf` and `_apply_quality_filter` in [variant_annotation.py](../../src/bioagent/pipelines/variant_annotation.py). Detection is in [detector.py](../../src/bioagent/agent/detector.py) (`##fileformat=VCF`).

**Likely interview questions**
1. *What are the first five columns of a VCF?* — CHROM, POS, ID, REF, ALT. Then QUAL, FILTER, INFO.
2. *What does QUAL mean?* — A Phred-scaled confidence that a variant exists at that site: −10·log10(probability the call is wrong).
