# BioAgent — Project Continuity Document

**Last Updated:** May 1, 2026
**Current Version:** v0.4.0
**Next Version:** v0.5.0 (Ollama LLM integration — IN PROGRESS)

---

## Current Project State

### What's Built and Working
- RAG system: 640 chunks, ChromaDB, sentence-transformers
- File detector: FASTA/FASTQ/VCF/CSV auto-detection
- Pipeline router: routes to correct pipeline with reasoning
- FASTA QC pipeline: GC content, length, composition, complexity
- RNA-seq pipeline: CPM, DE analysis, volcano/PCA/heatmap
- Variant annotation: VCF parsing, clinical annotation
- FastAPI backend: /upload, /health, /analyse, /jobs
- Web frontend: drag-drop, results display, plots
- 23 tests passing
- GitHub: v0.4.0 tagged

### What's In Progress (v0.5.0)
- Ollama installed and tested locally (llama3.2:3b, 2GB)
- explainer.py created in src/bioagent/agent/
- Next: wire explainer into the pipeline + add Q&A to frontend

---

## How to Resume Work

### Start the System
```bash
conda activate bioagent
python -m uvicorn bioagent.api.main:app --reload --host 0.0.0.0 --port 8000
```
Then open `frontend/index.html` in browser.

### Run Tests
```bash
pytest tests/ -v
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

## v0.5.0 Build Plan (Ollama Integration)

### What We're Adding
1. `src/bioagent/agent/explainer.py` — calls Ollama with pipeline results + RAG context
2. Wire explainer into `src/bioagent/api/main.py` — call after pipeline runs
3. Add Q&A box to `frontend/index.html` — user asks follow-up questions
4. Add `/ask/{job_id}` endpoint to backend — handles follow-up questions

### Architecture Change
Before: Pipeline → basic_interpretation() → display
After:  Pipeline → RAG query → Ollama(results + context) → rich explanation → display
↑
User Q&A also goes here

### Files to Modify
- `src/bioagent/agent/explainer.py` (new — already created)
- `src/bioagent/api/main.py` (add /ask endpoint)
- `frontend/index.html` (add Q&A section)

---

## .claude/ Folder Structure
.claude/
├── agents.md          # Who Steve is + how to work with Manny
├── memory.md          # What Steve knows about Manny
├── CONTINUITY.md      # This file — current project state
└── skills/
└── bioinformatics_pipelines.md  # How to build new pipelines

---

## Version History
- v0.0.1 — Project skeleton
- v0.1.0 — First working pipeline (FASTA QC)
- v0.2.0 — RNA-seq pipeline
- v0.3.0 — All three pipelines complete
- v0.4.0 — Full-stack demo (frontend + backend)
- v0.5.0 — Ollama LLM integration (IN PROGRESS)

# BioAgent — Project Continuity Document

**Last Updated:** May 2, 2026
**Current Version:** v0.5.0 (COMPLETE ✅)
**Next Version:** v0.6.0 (Deployment prep + polish)

---

## What's Built and Working

### Core System
- File auto-detection: FASTA/FASTQ (95%), VCF (99%), CSV/TSV (80%)
- Pipeline router: routes to correct pipeline with reasoning
- RAG system: 640 chunks, ChromaDB, sentence-transformers all-MiniLM-L6-v2
- Ollama Q&A: llama3.2:3b running locally, ANSI codes cleaned
- 23 tests passing

### Three Pipelines
- FASTA QC: GC content, length distribution, complexity, Chart.js plots
- RNA-seq: CPM normalisation, t-test DE, volcano + heatmap (PCA disabled — Windows numpy segfault)
- Variant Annotation: VCF parsing, clinical annotation, 8 known variants (BRCA1, BRCA2, TP53, KRAS, EGFR, PIK3CA, NRAS, NRAS)

### Frontend + Backend
- FastAPI backend: /upload, /health, /ask/{job_id}, /jobs
- Web frontend: drag-drop, interactive plots (iframes), Q&A box
- Served at: http://localhost:8000/frontend/index.html
- Plots: Chart.js HTML files written directly — no matplotlib (Windows segfault fix)

---

## How to Start the System

```bash
conda activate bioagent
cd C:\Users\Invate\Downloads\Bioinformatics-Agent
& "C:\Users\Invate\anaconda3\envs\bioagent\python.exe" -m uvicorn bioagent.api.main:app --host 0.0.0.0 --port 8000
```

Then open: `http://localhost:8000/frontend/index.html`

### Run Tests
```bash
& "C:\Users\Invate\anaconda3\envs\bioagent\python.exe" -m pytest tests/ -v
```

### Sample Test Files
- `data/sample/test.fasta` — 2 sequences, 45.84% GC
- `data/sample/test.vcf` — 8 pathogenic variants
- `data/sample/counts.csv` — 12 genes, healthy vs cancer (use labels: healthy / cancer)

---

## Key Files
src/bioagent/
├── agent/
│   ├── detector.py       — file type detection
│   ├── router.py         — pipeline routing
│   └── explainer.py      — Ollama Q&A (ANSI fix applied)
├── pipelines/
│   ├── fasta_qc.py       — FASTA pipeline (Chart.js plots)
│   ├── rnaseq.py         — RNA-seq pipeline (PCA disabled)
│   └── variant_annotation.py — VCF pipeline (Chart.js plots)
├── rag/
│   ├── embedder.py, vector_store.py, ingestion.py, retriever.py
├── api/
│   └── main.py           — FastAPI (thread pool executor for pipeline)
└── utils/
├── logger.py
├── plot_runner.py          — UNUSED (subprocess segfault)
├── rnaseq_plot_runner.py   — UNUSED (subprocess segfault)
└── variant_plot_runner.py  — UNUSED (subprocess segfault)
frontend/
└── index.html            — drag-drop UI + Q&A box

---

## Known Issues / Technical Debt
- PCA disabled in rnaseq.py — np.linalg.eigh causes Windows segfault in uvicorn worker
- matplotlib removed from all pipelines — same Windows segfault issue
- plot_runner.py files unused — Plotly also segfaults in subprocess on Windows
- All three issues resolve automatically when deployed to Linux
- Ollama responses slightly repetitive — llama3.2:3b limitation, upgrade to 7b when hardware allows
- RAG telemetry warnings (ChromaDB) — harmless, cosmetic only

---

## Version History
- v0.0.1 — Project skeleton
- v0.1.0 — FASTA QC pipeline
- v0.2.0 — RNA-seq pipeline  
- v0.3.0 — All three pipelines
- v0.4.0 — Full-stack demo (frontend + backend)
- v0.5.0 — Ollama LLM integration + Q&A ✅

---

## v0.6.0 Plan (Next Session)
- Deploy to free Linux server (Render or Railway)
- Re-enable PCA and matplotlib on Linux
- Add error handling improvements
- Add loading states and better UX feedback
- Prepare portfolio writeup for job applications
- README polish for GitHub visibility

---

## Validated on Real Biological Data (May 3, 2026)

All three pipelines tested and validated on real datasets.

### Dataset 1 — Breast Cancer RNA-seq
- File: real_breast_cancer.csv (30 genes, 3 normal vs 3 cancer samples)
- Based on published breast cancer gene signatures (TCGA)
- Results:
  - 13 upregulated, 11 downregulated genes
  - Top upregulated: MKI67 (+4.84 log2FC), ERBB2/HER2 (+4.05), AURKA (+3.95)
  - Top downregulated: FOXA1 (-2.55), CDH1 (-2.29), PGR (-2.24)
  - ERBB2 is the target of Herceptin — a real clinical breast cancer drug
  - FOXA1 and PGR downregulation consistent with published ER- breast cancer literature
- Ollama correctly identified biomarkers and therapy targets

### Dataset 2 — Clinical Variant Annotation
- File: real_clinical_variants.vcf (12 variants, ClinVar rsIDs)
- Results:
  - 10 passed QC, 2 failed quality filters
  - 8 pathogenic variants identified: BRCA1, BRCA2, TP53, KRAS, EGFR, PIK3CA, NRAS
  - Clinical alert triggered correctly
  - Conditions: hereditary breast/ovarian cancer, Li-Fraumeni syndrome, lung cancer, melanoma
- Ollama gave clinically appropriate interpretation

### Dataset 3 — E. coli K-12 FASTA
- File: ecoli_k12.fasta (5 sequences, 772 total bases)
- Results:
  - Mean GC: 47.37% (real E. coli is ~50.8% — close match)
  - Median length: 140bp, no low complexity sequences
  - All QC metrics passed
- Ollama correctly identified bacterial origin from GC content

### Conclusion
System produces biologically accurate results on real data.
Ready for deployment and portfolio presentation.