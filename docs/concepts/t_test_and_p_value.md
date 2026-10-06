# t-test and p-value

**What it is:** A test of whether two groups' averages differ by more than you'd expect from random sampling noise alone. The p-value is the output.

**Biology / intuition:** Three healthy and three tumour samples will never give identical numbers, even for a gene that does nothing. The t-test asks: "Is the gap between group means big compared with the scatter inside each group?" A big gap with tight scatter gives a small p-value.

**Formula in words:** t = (difference in group means) ÷ (standard error of that difference). The p-value is the probability of seeing a t at least this extreme **if there were truly no difference** (the null hypothesis).

**What p is NOT:** It is not the probability the gene is truly changed, and not the probability you're wrong. p = 0.04 means "if nothing were going on, data this extreme would occur about 4% of the time".

**When it's wrong to use:**
- RNA-seq counts are not normally distributed. They are skewed and their variance rises with the mean, which is why DESeq2/edgeR use a negative binomial model.
- With n = 3 per group the variance estimate is very unstable. DESeq2 borrows information across genes to fix this.
- Student's t-test (scipy's default, `equal_var=True`) assumes both groups have equal variance. Welch's t-test (`equal_var=False`) does not and is the safer default.
- Run on thousands of genes, 5% of nulls will give p<0.05 by chance (see [multiple_testing_fdr.md](multiple_testing_fdr.md)).

**Where in our code:** `stats.ttest_ind(control_counts, treatment_counts)` in `_differential_expression`, [rnaseq.py](../../src/bioagent/pipelines/rnaseq.py). It is Student's, run on raw counts.

**Likely interview questions**
1. *What does p = 0.03 mean?* — If the gene had no real difference, data this extreme would appear about 3% of the time. It isn't a 97% chance the gene is real.
2. *Why isn't a t-test ideal for RNA-seq?* — Counts are discrete, skewed and mean-dependent, n is tiny, and I'm testing thousands of genes. A negative binomial model with shared dispersion estimation (DESeq2) fits the data better.
