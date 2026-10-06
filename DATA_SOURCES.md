# DATA_SOURCES — real public test data

Rule: validation (see the VALIDATION RULE in [CLAUDE.md](CLAUDE.md)) needs data a stranger can download. Record the exact accession, version and date downloaded. Never commit large downloads: put them under `data/raw/` (already git-ignored) and record how to re-fetch them here.

I have not downloaded any of these yet. Accessions and sizes below are from memory. **Verify each one against the source page before relying on it.**

| Data | Used for | Where to get it | Size / licence notes |
|---|---|---|---|
| **airway RNA-seq** (Himes et al. 2014; human airway smooth muscle cells, dexamethasone vs untreated, 4 cell lines × 2) | Validating RNA-seq DE against pydeseq2/DESeq2 (L2). Also the classic DESeq2 tutorial dataset, so published results exist to compare | Bioconductor `airway` package (R), or a counts CSV exported from it. GEO series GSE52778 | Small (about 64k genes × 8 samples). Public. Cite Himes et al. 2014, PLoS ONE |
| **Insulin mRNA NM_000207.3** | Biopython-validated FASTA tests (L3); ORF/CDS finder (CDS annotated at 60–392, 110-aa preproinsulin) | NCBI Nucleotide: search NM_000207.3, download FASTA and GenBank | Tiny. Public domain NCBI data |
| **ClinVar subset** | Validating variant annotation (rsID to gene and significance) | NCBI ClinVar FTP: `variant_summary.txt.gz` or the VCF. Filter to the genes of interest | Full file is large (hundreds of MB). Commit only a small filtered extract. NCBI data is public |
| **GIAB (Genome in a Bottle) VCF** (e.g. HG002) | A real, truth-set-backed VCF to test parsing at realistic scale | NIST GIAB FTP | Whole-genome VCFs are huge. Use a single chromosome or region slice. Public |
| **E. coli K-12 MG1655 genome** | Replace the 1 KB excerpt in `data/sample/ecoli_k12.fasta`; GC check against the known genome | NCBI RefSeq NC_000913.3. A compressed copy already sits in the repo as `data/sample/ecoli_k12.fna.gz` (about 1.4 MB, tracked) | About 4.6 Mb. Public. NC_000913.3 GC is roughly 50.8%, so verify by computing it with Biopython |

## Existing sample files, provenance status

| File | What it is | Status |
|---|---|---|
| `data/sample/counts.csv` | 12 genes, 3 vs 3 | Synthetic/toy. Fine for unit tests, not for claims |
| `data/sample/real_breast_cancer.csv` | 30 genes, 3 vs 3 | **Provenance unknown**, looks hand-assembled. Don't call it real public data until sourced |
| `data/sample/real_clinical_variants.vcf` | 12 variants | **Provenance unknown.** Annotation is checked only against our own 8-entry dict |
| `data/sample/ecoli_k12.fasta` | About 1 KB excerpt | An excerpt, so its GC% can't be compared to the whole genome |
| `data/knowledge/*.txt` | Your MSc lecture notes (RAG source) | Tracked in a public repo. Check your course's rules on sharing notes |
