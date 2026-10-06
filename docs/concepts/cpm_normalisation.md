# CPM normalisation (Counts Per Million)

**What it is:** Rescale each sample's gene counts so every sample is "per million reads", making samples comparable.

**Biology / intuition:** Sequencers don't read the same number of fragments from every sample. A sample sequenced to 20M reads will show roughly double the counts of one sequenced to 10M, even if the biology is identical. CPM removes that depth difference.

**Formula in words:** for each gene in a sample, divide its count by the sample's total counts, then multiply by 1,000,000.

**When it's wrong to use:**
- It doesn't correct for gene length, so you can't compare gene A to gene B within a sample (that needs TPM/FPKM).
- It assumes most genes don't change. If a few huge genes swing (or a sample is dominated by one transcript), every other gene looks shifted. DESeq2's median-of-ratios and edgeR's TMM are more robust to this ("composition bias").
- Don't feed CPM to count-based models like DESeq2. They want raw counts and do their own normalisation.

**Where in our code:** `_normalise_cpm` in [rnaseq.py](../../src/bioagent/pipelines/rnaseq.py). **Known issue:** the CPM table is only used for the heatmap. The t-test uses raw counts (see [DECISIONS.md](../../DECISIONS.md)).

**Likely interview questions**
1. *Why normalise at all?* — Different sequencing depth per sample would make a gene look up- or down-regulated when it isn't. CPM scales every sample to the same depth.
2. *Is CPM enough for differential expression?* — No. It fixes depth, not composition bias, and it isn't a statistical model. Production tools (DESeq2, edgeR) use their own size factors and a negative binomial model. I'm validating my pipeline against pydeseq2 for that reason. *(Earn in W3/L2.)*
