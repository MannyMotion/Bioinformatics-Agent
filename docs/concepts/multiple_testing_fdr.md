# Multiple testing and FDR (Benjamini-Hochberg)

**What it is:** When you run thousands of tests, many "significant" results are flukes. FDR control limits the *share* of your significant calls that are false.

**Biology / intuition:** Test 20,000 genes at p<0.05 and, even if nothing changed at all, you'd expect about 1,000 false hits. A volcano plot full of "significant" genes can be mostly noise. FDR asks: "Of the genes I call significant, what fraction do I expect to be false?"

**How Benjamini-Hochberg works in words:**
1. Sort all p-values smallest to largest (rank 1…m).
2. For each p-value, compute p × m ÷ rank. This is its adjusted value.
3. Make sure adjusted values never decrease as you go to bigger p (take a running minimum from the bottom).
4. Genes with adjusted p (called q or padj) below 0.05 are your significant set. At 5% FDR, about 5% of that set is expected to be false.

**When it's wrong / limits:**
- It controls the *expected proportion* of false discoveries, not the chance of any single one (that is Bonferroni, which is stricter).
- It assumes tests are independent or positively dependent. Genes are correlated, but BH is generally fine in practice.
- It can't rescue a bad underlying test. If the p-values come from the wrong model, BH adjusts wrong numbers.

**Where in our code:** **Not implemented yet.** `significant` uses raw p<0.05. Lesson L1 adds it. `statsmodels.stats.multitest.multipletests(p, method="fdr_bh")` is the reference to test against.

**Likely interview questions**
1. *Why not just use p<0.05 for each gene?* — With 20,000 tests I'd expect about 1,000 false positives by chance. I correct with BH so the false discovery rate among my hits is controlled at 5%.
2. *Difference between Bonferroni and BH?* — Bonferroni controls the chance of even one false positive and is very strict. BH controls the proportion of false positives among calls, so it keeps more power. That's the usual choice in genomics.
