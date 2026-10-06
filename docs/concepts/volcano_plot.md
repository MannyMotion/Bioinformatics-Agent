# Volcano plot

**What it is:** A scatter plot where every dot is a gene. x = log2 fold change (how big the change is), y = −log10(p-value) (how confident the test is).

**Biology / intuition:** You want genes that changed *a lot* and where the test says it's *unlikely to be chance*. Those sit in the upper left (down-regulated) and upper right (up-regulated) corners. The "volcano" shape appears because big fold changes tend to come with small p-values, and nothing happens in the middle (low change, low significance).

**Why −log10(p):** p = 0.05 gives 1.3, p = 0.001 gives 3, p = 1e-10 gives 10. Small p-values become tall points, so "more significant" means "higher up".

**How to read it:** Vertical lines at ±1 (2× change) and a horizontal line at the p-threshold. Coloured points pass both cuts.

**When it's misleading:**
- Without FDR correction, the y-axis overstates significance across thousands of tests.
- Low-count genes can have huge fold changes and still look interesting.
- With tiny n, a "perfect" volcano may reflect unstable variance estimates.

**Where in our code:** `_generate_plots` in [rnaseq.py](../../src/bioagent/pipelines/rnaseq.py), Chart.js scatter. Note it plots `-log10(p + 1e-10)`, so p = 0 is capped at about 10.

**Likely interview questions**
1. *What does a gene in the top right mean?* — Strongly up-regulated with a small p-value: a candidate, which I'd confirm with FDR and independent data.
2. *Why −log10?* — It flips and stretches p-values so tiny p-values become tall, readable points.
