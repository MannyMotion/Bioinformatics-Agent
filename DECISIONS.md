# DECISIONS — BioAgent design decision log

One entry per decision: date, decision, why, alternatives, and **how I'd explain it in an interview**.
Entries marked *reconstructed* were written after the fact from the code and git history. The "why" is our best reading, not a recorded quote. Check each one with Manny as he works through the matching lesson.
New decisions go at the **top** of the log. Big or hard-to-reverse ones go through the Advisor skill first.

---

## 2026-10-06 — Teaching mode and the validation rule
- **Decision:** From now on Manny types all core (bioinformatics/stats/SQL/test) code. Claude teaches and reviews. No pipeline is "done" until it matches a trusted reference tool on real public data, with a test.
- **Why:** Much of the code was AI-written and Manny can't yet defend it in an interview. A portfolio is worth only what its owner can explain. Proof of correctness is what impressed people on the BioSeq Analyzer post (see [IDEAS.md](IDEAS.md)).
- **Alternatives:** Keep having Claude write code (faster, but doesn't build the skill or the story). Add more features first (less credible than fewer, validated ones).
- **Interview wording:** "I realised I'd moved faster than my understanding, so I rebuilt my workflow: I rewrite and validate each method myself against a reference tool, and I keep a log of what I can and can't explain."

## 2026-10-06 — Finding: the t-test runs on raw counts, not CPM
- **Decision/status:** Recorded as a known defect, not yet fixed. `_differential_expression` receives `filtered_df` (raw counts). CPM is only used for the heatmap. Earlier notes saying "t-test on CPM" were wrong.
- **Why it matters:** Different library sizes distort raw-count comparisons. This is exactly what normalisation is for, and it is currently bypassed for the statistics.
- **Next:** Lessons W3, L1, L2 address it. Replace with a validated method (pydeseq2) and compare.
- **Interview wording:** "My first version normalised for the heatmap but tested raw counts. I found it by reading my own code, and fixed it by validating against pydeseq2."

---

## Reconstructed decisions (from code + git history)

### ~2026-04 — FastAPI backend + vanilla JS frontend *(reconstructed)*
- **Decision:** FastAPI + Uvicorn for the API, a single HTML/CSS/vanilla-JS page for the UI.
- **Why:** FastAPI is lightweight, async, typed, and auto-documents. Vanilla JS avoids a build step and a framework while the focus is the bioinformatics.
- **Alternatives:** Streamlit (faster, less control, how the BioSeq Analyzer was built), Flask, React.
- **Interview wording:** "I wanted a real API that other tools could call, so FastAPI. I kept the front end deliberately simple because the value is in the analysis, not the UI framework."

### ~2026-05 — Local Ollama LLM (llama3.2:3b) *(reconstructed)*
- **Decision:** Explanations and Q&A come from a local model via the `ollama` CLI.
- **Why:** Free, no API key, and data never leaves the machine, which matters for patient or genomic data. Fits a 16 GB laptop.
- **Alternatives:** A hosted API (better quality, costs money, data leaves the machine).
- **Known downside:** The 3B model is weak. It also can't run on free hosting, so a public live demo needs pre-written explanations (Lesson L4). The model only receives summary counts, not gene lists, so it can't ground statements in specific results.
- **Interview wording:** "Privacy and cost drove it. Genomic data is sensitive, so I kept inference local, and I'm honest that a small model is a trade-off."

### ~2026-04 — ChromaDB RAG over my MSc notes *(reconstructed)*
- **Decision:** Chunk my MSc lecture notes (500 chars, 100 overlap), embed with all-MiniLM-L6-v2, store in ChromaDB, retrieve the top 2–3 chunks into the LLM prompt.
- **Why:** Gives the LLM domain context from material I actually studied, and demonstrates RAG, a skill employers ask about.
- **Alternatives:** Fine-tuning (expensive), no retrieval (the model guesses).
- **Known caveat:** The README says "grounded, not hallucinated". That is too strong. Retrieval reduces but doesn't eliminate hallucination.
- **Interview wording:** "RAG lets a small local model answer using my own course notes instead of only its training data. It reduces hallucination. It doesn't prove correctness."

### ~2026-07 — SQLite job store *(reconstructed)*
- **Decision:** Persist each job's JSON in `jobs.db` (table `jobs`: job_id, data, created_at). A fresh connection is opened per call.
- **Why:** The earlier in-memory dict lost jobs on restart and grew without bound. SQLite is in the standard library: no extra service, no cost.
- **Alternatives:** In-memory dict, Postgres (overkill now), files on disk.
- **Interview wording:** "Zero-dependency persistence so follow-up Q&A survives restarts. I'd move to Postgres when there are multiple users."

### ~2026-07 — Two-step agentic flow: analyst.py then pipelines *(reconstructed)*
- **Decision:** Replace the fixed router with a two-phase flow. `/upload` detects and profiles the file, infers groups and proposes questions. `/analyze` runs the chosen pipeline.
- **Why:** Lets users start without knowing which tool to run, and lets the system infer control vs treatment from column names instead of fixed prefixes.
- **Honest note:** The "agent" is deterministic rules (keyword matching on column names), not an LLM planner. That makes it fast, offline and predictable, but I shouldn't oversell it.
- **Interview wording:** "The profiling step is rule-based so it's instant and testable. The LLM only explains results. I'd say it's agent-like workflow, not an autonomous agent."

### ~2026-04 — Differential expression: two-sample t-test *(reconstructed, corrected)*
- **Decision:** Per-gene `scipy.stats.ttest_ind` between groups, with log2FC from group means (+1 pseudo-count). Significant if p<0.05 and |log2FC|≥1.
- **Why (as recorded in code):** "dependency simplicity while learning the concepts." DESeq2 is acknowledged as the production standard.
- **Corrections:** It runs on **raw counts**, not CPM. It is Student's (equal-variance) t-test. There is no multiple-testing correction.
- **Alternatives:** DESeq2/edgeR (negative binomial), pydeseq2 (Python), limma-voom.
- **Next:** Lessons L1 (BH-FDR) and L2 (validate vs pydeseq2).
- **Interview wording:** "It was a learning-stage placeholder. I know why it's inadequate: counts aren't normal, n is tiny, and I test thousands of genes. That's why I'm validating against pydeseq2."

### ~2026-04 — Chart.js HTML plots instead of matplotlib *(reconstructed)*
- **Decision:** Plots are standalone HTML files using a vendored Chart.js (`static/chart.umd.min.js`).
- **Why:** matplotlib (and numpy eigh for PCA) segfaulted inside the Windows uvicorn worker. Vendoring Chart.js keeps plots working offline.
- **Alternatives:** Plotly (tried first, crashed), running matplotlib in a subprocess (tried, left unused `*_plot_runner.py`).
- **Cost:** PCA is disabled. Plots are less publication-grade.
- **Interview wording:** "I hit a platform-specific crash, isolated it, and chose a browser-rendered library rather than fight it."

### ~2026-04 — Curated dictionary instead of live ClinVar *(reconstructed)*
- **Decision:** `KNOWN_VARIANTS` is a hardcoded dict of 8 rsIDs mapped to gene/condition/significance.
- **Why:** Simple, offline, enough to demo the annotation workflow.
- **Caveat:** It is not ClinVar. The mappings are unverified against ClinVar, so the validation rule applies (Lesson on ClinVar, see [ROADMAP.md](ROADMAP.md)).
- **Interview wording:** "It demonstrates the workflow. For real use you'd query ClinVar or run VEP. I'm validating my mappings against ClinVar."

### ~2026-04 — Content-based file detection *(reconstructed)*
- **Decision:** Read only the first 20 lines and decide type by structure (`@`+`+` FASTQ, `>` FASTA, `##fileformat=VCF`, delimiter consistency).
- **Why:** Extensions are unreliable, and FASTQ files can be 50 GB.
- **Interview wording:** "I never read a whole file for detection, and I trust the content over the extension."
