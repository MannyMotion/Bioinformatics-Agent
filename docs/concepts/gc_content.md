# GC content

**What it is:** The percentage of bases in a sequence that are G or C.

**Biology / intuition:** G–C pairs have three hydrogen bonds, A–T have two, so GC-rich DNA is more stable. GC varies a lot between species (a few % in some AT-rich parasites up to over 70% in some bacteria; human is roughly 41%; *E. coli* K-12 about 50.8%). Within a genome it varies too: genes and CpG islands are richer than gaps. Unexpected GC in a sample can signal contamination (another organism mixed in) or sequencing bias.

**Formula in words:** (number of G + number of C) ÷ (number of A + C + G + T) × 100. We leave out N (unknown bases) from the denominator.

**When it's wrong to use:**
- A single GC% for a whole file hides variation: one contaminant contig can be invisible in the average.
- The mean of per-sequence GC% (what we compute) treats a 100 bp and a 1 Mb sequence equally. It is **not** the genome GC. Weight by length for that.
- Comparing a short excerpt to a published genome-wide value (as the README did for *E. coli*) is invalid.
- Our thresholds (35% / 65%) are generic QC flags, not species-specific truths.

**Where in our code:** `_analyse_sequence` and `_generate_warnings` in [fasta_qc.py](../../src/bioagent/pipelines/fasta_qc.py). Validate against Biopython's `gc_fraction` (Lesson L3).

**Likely interview questions**
1. *Why does GC content matter in QC?* — It shows species/contamination clues and affects sequencing coverage and PCR/amplification bias. A sample whose GC doesn't match the expected organism deserves a closer look.
2. *Is your average GC the genome's GC?* — Not exactly. I average per-sequence values, so I'd weight by length to get the true genome figure. I caught that while validating.
