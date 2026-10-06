# INTERVIEW BANK — BioAgent

Answers are in Manny's voice. Tag meaning: **[earned]** = Manny has typed and explained this, so the answer is safe to say. **[learn in X]** = not yet earned. Don't say it in an interview until the lesson is done and the LEARNING_LOG row says "yes". Many answers are deliberately honest about weaknesses. Admitting a weakness and showing how you'd fix it reads as seniority.

Nothing is earned yet (teaching mode starts at W1).

## 1. The 30-second pitch

**Q: Tell me about BioAgent.**
"It's a tool where a biologist uploads a FASTA, a gene-count table or a VCF, and it detects the file type, runs the matching analysis (sequence QC, differential expression, or variant annotation), and explains the result in plain English using a local LLM grounded in my MSc notes. I built it to practise end-to-end bioinformatics engineering, and I'm now validating each pipeline against reference tools like pydeseq2 and Biopython." **[learn in W1–W5; the validation claim becomes true after L1–L3]**

## 2. Design

**Q: Why FastAPI?** "A real API that other tools can call, typed and auto-documented. I kept the front end simple so the effort went into the analysis." **[learn in W5]**

**Q: Is it really an 'agent'?** "I'd call it an agent-like workflow. A deterministic analyst step profiles the data, infers groups from column names and proposes questions. The LLM only explains the results. That's deliberate: the analysis stays testable and offline." **[learn in W5]**

**Q: Why SQLite for jobs?** "Zero-dependency persistence so follow-up questions survive a restart. If there were many users I'd move to Postgres." **[learn in W5]**

**Q: Why a local LLM?** "Privacy and cost. Genomic data is sensitive, so I keep inference local. I accept lower quality for that." **[learn in W5]**

## 3. Biology and statistics

**Q: Walk me through your RNA-seq pipeline.** "Load the count matrix, drop very low-count genes, normalise for library size, compare two groups per gene, compute log2 fold change and p-values, correct for multiple testing, and plot a volcano and heatmap." *(True once L1 is done. Today the pipeline lacks the correction.)* **[learn in W3, L1]**

**Q: What's CPM and why use it?** "Counts per million: each count divided by the sample's total, times a million. It removes sequencing-depth differences so samples are comparable." See [cpm_normalisation.md](../concepts/cpm_normalisation.md). **[learn in W3]**

**Q: What does log2 fold change mean?** "log2 of treatment over control. +1 means doubled, −1 means halved, +4 means 16×." See [log2_fold_change.md](../concepts/log2_fold_change.md). **[learn in W3]**

**Q: Explain a p-value.** "If there were truly no difference, how often would I see data this extreme? It's not the probability that the gene is real." See [t_test_and_p_value.md](../concepts/t_test_and_p_value.md). **[learn in W3]**

**Q: Why not just threshold p<0.05?** "With thousands of genes you get hundreds of false positives by chance. I use Benjamini-Hochberg to control the false discovery rate." See [multiple_testing_fdr.md](../concepts/multiple_testing_fdr.md). **[learn in L1]**

**Q: Is a t-test right for RNA-seq?** "Not really. Counts are skewed, their variance rises with the mean, and n is tiny. DESeq2 uses a negative binomial model and shares information across genes. My first version used a t-test as a learning placeholder, and I validated it against pydeseq2." **[learn in W3, L2]**

**Q: Your t-test used raw counts, not CPM. Why?** "That was a flaw. I normalised for the heatmap but tested raw counts. I found it reading my own code and it's logged in my decisions file. The validated version uses pydeseq2's own normalisation." **[learn in W3, L2]**

**Q: What's GC content and what does it tell you?** "The fraction of G and C bases. It varies by species and can reveal contamination or bias." See [gc_content.md](../concepts/gc_content.md). **[learn in W2]**

**Q: What do the first columns of a VCF mean? What's QUAL?** "CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO. QUAL is a Phred-scaled confidence that the variant exists." See [vcf_format.md](../concepts/vcf_format.md). **[learn in W4]**

**Q: How do you decide if a variant is pathogenic?** "I don't decide: I look it up. Classification follows ACMG/AMP guidelines, and ClinVar aggregates submitters' classifications. My pipeline currently uses a small curated table and I'm validating it against a real ClinVar extract." **[learn in W4]**

## 4. Honesty and "what would you improve"

**Q: What's the weakest part?** "The statistics were a t-test with no multiple-testing correction, my annotation was a hardcoded 8-variant table, and my sample data was small and toy-like. I've logged each, and I'm fixing them in order of credibility: FDR, then validation against pydeseq2, then Biopython and ClinVar checks." **[earned when L1–L3 are done]**

**Q: Your low-complexity check flags every sequence above ~860 bp. Did you know?** "Yes: only 256 possible 4-mers, so the unique-to-total ratio is capped. I'd switch to a windowed method like DUST." **[learn in W2]**

**Q: How do you know your output is correct?** "I compare it to a trusted tool on public data and have a test that fails if they disagree: pydeseq2 for RNA-seq, Biopython for sequences, ClinVar for variants." **[learn in L2, L3]**

**Q: What would you do with only 3 samples per group?** "Treat results as hypotheses, not findings. With n = 3 variance is unstable. I'd say so, use a method that shares dispersion information (DESeq2), and plan for more replicates." **[learn in W3]**

## 5. Business

**Q: How could a company use this?** "A small lab without a bioinformatician could get a first-pass QC and a sanity-checked DE list in minutes, with plain-English notes. For a core facility it could triage submissions (is the FASTA contaminated? did the RNA-seq groups separate?) before a specialist spends time." Caveats: needs validation against their pipelines, audit trails, and data-governance review. **[earned after L2 + L5]**

**Q: How would you deploy it?** "Pipelines on a small cloud service, local-LLM explanations replaced by pre-written or hosted-model text, and auth plus rate limits. Containerised with Docker." **[learn in L4]**
