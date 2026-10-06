# ROADMAP — lessons in order

**Current lesson: W1** (update this line when a lesson finishes).
Each lesson follows [.claude/skills/tutor/SKILL.md](.claude/skills/tutor/SKILL.md). Each ends with Manny explaining it back. Mastery over pace: 2 weeks behind with real understanding beats being on time without it.

Status key: ☐ not started · ◐ in progress · ☑ done and explained back

## PHASE 0 — Understand what you have (walkthroughs, no new features)

| # | Lesson | You will be able to explain | Status |
|---|---|---|---|
| W1 | File detection ([detector.py](src/bioagent/agent/detector.py)) | How a file is recognised from its first lines; FASTA vs FASTQ vs VCF structure; Phred encoding | ☐ |
| W2 | FASTA QC and GC-content biology ([fasta_qc.py](src/bioagent/pipelines/fasta_qc.py)) | What GC% is and why it differs between species; low-complexity; why an unweighted mean of per-sequence GC is not genome GC | ☐ |
| W3 | RNA-seq DE ([rnaseq.py](src/bioagent/pipelines/rnaseq.py)) | CPM, log2FC, t-test, p-values, why the volcano plot looks like that, and the raw-counts-vs-CPM bug | ☐ |
| W4 | Variant annotation ([variant_annotation.py](src/bioagent/pipelines/variant_annotation.py)) | VCF fields, ClinVar, what "pathogenic" means, and why a hardcoded dict is not ClinVar | ☐ |
| W5 | The agent layer ([analyst.py](src/bioagent/agent/analyst.py), RAG, Ollama, [job_store.py](src/bioagent/api/job_store.py)) | Request flow end to end; what is rules and what is LLM; first SQL queries on jobs.db | ☐ |

## PHASE 1 — Credibility (the validation rule in action)

| # | Lesson | Done when | Status |
|---|---|---|---|
| L1 | Benjamini-Hochberg FDR on the t-test | You type the BH correction and a test checks it against `statsmodels.multipletests` | ☐ |
| L2 | Validate RNA-seq DE against pydeseq2 on a public dataset (e.g. airway) | A test shows the gene lists agree within a stated tolerance | ☐ |
| L3 | Biopython-validated tests for the FASTA pipeline | GC, length and counts match Biopython on a real NCBI sequence | ☐ |

(Variant validation against ClinVar is queued in [IDEAS.md](IDEAS.md) and becomes a lesson after L3.)

## PHASE 2 — Show it

| # | Lesson | Done when | Status |
|---|---|---|---|
| L4 | Free live demo: pipelines online plus pre-written explanations for the example files (Ollama can't run on free hosting) | A public URL works | ☐ |
| L5 | README with one headline validated result | A recruiter can read the proof in 30 seconds | ☐ |
| L6 | LinkedIn post | Posted, leading with proof, not features | ☐ |

## PHASE 3 — Grow it

Pick from [IDEAS.md](IDEAS.md), highest hiring value and learning value first.
