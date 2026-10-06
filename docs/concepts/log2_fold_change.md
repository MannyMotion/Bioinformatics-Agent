# log2 fold change (log2FC)

**What it is:** How much a gene's expression changed between two groups, on a log base-2 scale.

**Biology / intuition:** "ERBB2 is 16× higher in tumour" is the same as log2FC = +4. Doubling is +1, halving is −1, no change is 0. The log scale makes up and down symmetric (2× up = +1, 2× down = −1), which a plain ratio does not (2 vs 0.5).

**Formula in words:** log2 of (mean treatment expression ÷ mean control expression). We add 1 to each mean first (a pseudo-count) so we never take log of zero.

**When it's wrong to use:**
- Low-count genes have noisy, inflated fold changes (1 vs 3 counts is "log2FC 1.6" but meaningless). DESeq2 applies "shrinkage" to tame this.
- A big fold change with a tiny sample is a hypothesis, not a finding.
- The pseudo-count (+1) distorts genes with small counts.

**Where in our code:** `_differential_expression` in [rnaseq.py](../../src/bioagent/pipelines/rnaseq.py): `np.log2((mean_treatment + 1) / (mean_control + 1))`, computed on raw-count means. Threshold: |log2FC| ≥ 1 (2× change).

**Likely interview questions**
1. *What does log2FC = −2 mean?* — Expression in treatment is a quarter of control (2⁻² = 0.25).
2. *Why log2 and not the raw ratio?* — Symmetry between up and down, and fold changes behave additively on a log scale. A raw ratio is lopsided (up is 1 to infinity, down is 0 to 1).
