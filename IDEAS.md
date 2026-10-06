# IDEAS — backlog

Status: `idea` / `queued` / `in progress` / `done` / `dropped`. Effort: S (hours), M (days), L (week+).
**Lesson from Damilare Kareem's "BioSeq Analyzer" LinkedIn post** (Streamlit app, 15 tests validated against Biopython; an unannotated insulin mRNA NM_000207.3 gave a CDS at exactly NCBI's 60–392, translated to the 110-aa preproinsulin; ORFs in 6 frames; restriction mapping; commenters asked for intron/exon detection). What impressed people was **proof it's right**, not features. Rank everything below by that.

| # | Idea | What | Why (biology + hiring value) | Effort | Learning value for Manny | Status |
|---|---|---|---|---|---|---|
| 1 | DESeq2-grade differential expression | Replace the t-test with pydeseq2, keep a t-test+BH baseline for comparison | The #1 credibility gap. "Validated against pydeseq2" is a hireable sentence | M | Very high: stats, pandas, testing | queued (L1, L2) |
| 2 | Biopython-validated sequence core | Reimplement GC, length, composition and compare to Biopython on NM_000207.3 | Direct copy of what impressed people on the BioSeq post | S | High: loops, testing | queued (L3) |
| 3 | ORF finder / CDS finder, 6 frames | Find ORFs, translate, compare to NCBI's annotated CDS (insulin 60–392) | Central dogma made visible. Easy to prove right | M | High: loops, strings, codon table | idea |
| 4 | Real ClinVar-backed variant annotation | Replace the 8-rsID dict with a downloaded ClinVar subset (VCF/TSV) and join on rsID/position | The current annotation isn't ClinVar. This makes it true and testable | M | High: pandas merge, SQL | idea |
| 5 | Replace hand-made sample data with public data | airway (RNA-seq), GIAB or ClinVar subset (VCF), full E. coli genome | "Real data" is currently toy. Public data is citable | S | Medium | queued (see [DATA_SOURCES.md](DATA_SOURCES.md)) |
| 6 | Give the LLM the actual results | Pass top genes and p-values (not just counts) into the prompt, and show them in the answer | The LLM currently sees only counts, so it can't be grounded in the results | S | Medium | idea |
| 7 | PCA and sample QC | Re-enable PCA (Linux/hosting) and add sample-correlation QC | Reviewers expect PCA to check batch/outliers | M | Medium: linear algebra intuition | idea |
| 8 | Restriction-site mapping | Find enzyme cut sites in a sequence | Classic molecular-biology tool, easy to validate vs Biopython `Restriction` | S | Medium | idea |
| 9 | Intron/exon awareness | Handle GenBank features so exon/CDS boundaries come from annotation | Asked for by commenters on the BioSeq post | M | Medium: parsing GenBank | idea |
| 10 | Pathway / GO enrichment | Test DE gene lists for enriched pathways (hypergeometric test, FDR) | Turns a gene list into biology. Core RNA-seq skill | M | High: stats | idea |
| 11 | Versioned results and a job history UI | Query the SQLite job store, list past runs | Shows SQL and product thinking | S | High: SQL | idea |
| 12 | Free public demo | Host pipelines online with pre-written explanations | Recruiters click links, not clone repos | M | Medium: deployment | queued (L4) |
| 13 | Proper test suite | Assert correct numbers (known GC of a known sequence, known p-value) not just "it runs" | The current tests are smoke tests | S | Very high | queued (part of L1–L3) |
| 14 | Clean the repo | Remove unused `*_plot_runner.py`, stray matplotlib imports, untrack outputs PNGs | Reads as polished to a reviewer | S | Low | idea (needs Manny's yes) |
