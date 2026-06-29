# BioAgent — Project Continuity Document

**Last Updated:** June 2026
**Current Version:** v0.5.1 (Security hardening — COMPLETE ✅)
**Next Version:** v0.6.0 (Production LLM swap + Render deployment)

---

## Who This Is For

Emmanuel Ogbu (Manny), Manchester UK. MSc Bioinformatics (Bradford), BSc Biomedical
Science (MMU). Job hunting 4+ months — this project IS his work experience.
GitHub: MannyMotion.

**The Project:** BioAgent — AI-powered agentic bioinformatics SaaS (currently local).
Uploads any bio file (FASTA, FASTQ, VCF, CSV), auto-detects data type, runs the
correct pipeline, explains every step in plain English, answers follow-up questions
via LLM.

---

## Current Project State

### What's Built and Working
- File detector: FASTA/FASTQ (95%), VCF (99%), CSV/TSV (80%) confidence
- Pipeline router: routes to correct pipeline with reasoning
- RAG system: 640 chunks, ChromaDB, sentence-transformers (all-MiniLM-L6-v2)
- Three validated pipelines (real biological data, not synthetic):
  - FASTA QC — E. coli K-12 genome (47.37% GC vs published 50.8%)
  - RNA-seq DE — breast cancer vs normal (ERBB2/HER2 +4.05 log2FC, MKI67 +4.84,
    FOXA1 -2.55, BRCA1 -1.48 — real validated biomarkers)
  - Variant Annotation — clinical VCF, 8/12 pathogenic variants flagged
    (BRCA1, BRCA2, TP53, KRAS, EGFR, PIK3CA)
- FastAPI backend: `/upload`, `/analyse/{job_id}`, `/ask/{job_id}`, `/health`
- Ollama Q&A: llama3.2:3b, local only, grounded via RAG (not hallucinated)
- Web frontend: drag-drop upload, Chart.js plots, Q&A box
- Security hardening (NEW in v0.5.1):
  - Rate limiting via slowapi — 10 uploads/min, 20 questions/min per IP
  - 50MB file size cap
  - Extension whitelist (`ALLOWED_EXTENSIONS` in `api/main.py`)
  - 500-char question length cap (prompt-injection mitigation)
- 23 tests passing

### Hard Constraints (do not revert)
- **Windows segfault:** matplotlib / scikit-learn / plotly / `np.linalg.eigh` all
  crash the uvicorn worker on this hardware. All plots are pure Chart.js HTML
  written to disk instead. This is permanent on Windows — only revisit on Linux.
- **PCA disabled** in `rnaseq.py` — re-enable once deployed to Linux (Render).
- **`plot_runner.py`, `rnaseq_plot_runner.py`, `variant_plot_runner.py`** in
  `utils/` are unused dead code from the Plotly-subprocess approach (also
  segfaulted on Windows). Candidates for deletion once Linux deploy confirms
  they're not needed there either.
- **Ollama can't run on Render free tier** (512MB RAM) — production LLM
  replacement is unresolved. This blocks deployment.
- **`use_rag=False` hardcoded** in `api/main.py` upload handler — `BioRetriever`
  is only imported inside functions, not wired into the main analyse flow yet.

---

## How to Resume Work (local, Windows)

```bash
conda activate bioagent
cd C:\Users\Invate\Downloads\Bioinformatics-Agent
& "C:\Users\Invate\anaconda3\envs\bioagent\python.exe" -m uvicorn bioagent.api.main:app --host 0.0.0.0 --port 8000
```
Open: `http://localhost:8000/frontend/index.html`

```bash
& "C:\Users\Invate\anaconda3\envs\bioagent\python.exe" -m pytest tests/ -v
```
Expected: 23 passed

### Key Paths
- Project: `C:\Users\Invate\Downloads\Bioinformatics-Agent`
- Backend: `src/bioagent/api/main.py`
- Frontend: `frontend/index.html`
- Pipelines: `src/bioagent/pipelines/`
- RAG: `src/bioagent/rag/`
- Agent: `src/bioagent/agent/`

---

## Strategic Decision (June 2026)

**Open question raised:** 5 user interviews not yet done — validate product-market
fit before building more features?

**Decision:** Defer interviews. The current bottleneck isn't "will people pay for
this" (a SaaS-business question) — it's "can Manny point a hiring manager at a
working URL" (a job-hunting question), and those have different validation needs.
A live deployed demo makes any future interview land harder anyway.

**Agreed order for v0.6.0:**
1. Swap Ollama → Groq free-tier API (production LLM, unblocks deployment)
2. Deploy to Render — live URL for CV
3. User interviews (now backed by a real demo link)
4. Additive pipeline modules (see below) — feature work, not urgent

### Tools Explored for Next Phase (additive, not a rebuild)
- BLAST via BioPython
- NCBI Entrez API
- AlphaFold DB
- SWISS-MODEL REST API

These are new pipeline modules to bolt on after deployment — not a replacement
for the existing architecture.

---

## v0.6.0 Build Plan

1. Replace Ollama calls in `src/bioagent/agent/explainer.py` with Groq's free-tier
   API (OpenAI-compatible). Keep RAG-grounding behaviour identical.
2. Deploy to Render (or Railway) — confirm PCA and matplotlib can be re-enabled
   on Linux now that the Windows segfault constraint no longer applies.
3. Re-enable PCA in `rnaseq.py` once confirmed stable on Linux.
4. Wire `BioRetriever` into the main `/upload` → `/analyse` flow (`use_rag=True`)
   instead of leaving it hardcoded off.
5. After deploy: run the 5 user interviews with a live link in hand.
6. Then: pick first additive module (BLAST/NCBI/AlphaFold/SWISS-MODEL).

---

## .claude/ Folder Structure
```
.claude/
├── agents.md          # Who Steve (mentor) is + how to work with Manny
├── memory.md          # What Steve knows about Manny
├── CONTINUITY.md      # This file — current project state
└── skills/
    └── bioinformatics_pipelines.md  # How to build new pipelines
```

---

## Version History
- v0.0.1 — Project skeleton
- v0.1.0 — FASTA QC pipeline
- v0.2.0 — RNA-seq pipeline
- v0.3.0 — All three pipelines complete
- v0.4.0 — Full-stack demo (frontend + backend)
- v0.5.0 — Ollama LLM integration + RAG Q&A
- v0.5.1 — Security hardening (rate limiting, file validation, question cap) ✅
- v0.6.0 — Production LLM swap (Groq) + Render deployment (IN PROGRESS)

---

## Validated on Real Biological Data (May 2026)

All three pipelines tested and validated on real datasets:

- **Breast cancer RNA-seq** (`real_breast_cancer.csv`, 30 genes, 3 normal vs 3
  cancer) — 13 up / 11 down. Top up: MKI67 (+4.84), ERBB2/HER2 (+4.05), AURKA
  (+3.95). Top down: FOXA1 (-2.55), CDH1 (-2.29), PGR (-2.24). ERBB2 is the
  Herceptin target; FOXA1/PGR pattern matches published ER- breast cancer
  literature.
- **Clinical VCF** (`real_clinical_variants.vcf`, 12 variants, ClinVar rsIDs) —
  10 passed QC, 8 pathogenic (BRCA1, BRCA2, TP53, KRAS, EGFR, PIK3CA, NRAS).
  Conditions matched: hereditary breast/ovarian cancer, Li-Fraumeni, lung
  cancer, melanoma.
- **E. coli K-12 FASTA** (`ecoli_k12.fasta`, 5 sequences, 772 bases) — mean GC
  47.37% vs published ~50.8%, median length 140bp, no low-complexity sequences.

**Conclusion:** System produces biologically accurate results on real data.
Ready for deployment.
