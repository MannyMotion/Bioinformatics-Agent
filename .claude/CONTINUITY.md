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