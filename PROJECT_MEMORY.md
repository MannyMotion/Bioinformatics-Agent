# 🧬 BioAgent — Living Project Memory

> **This is the single source of truth for where the project is.**
> Read it at the start of every session. Update it at the end of every session.
> It supersedes the older `.claude/CONTINUITY.md` (which is now historical).

**Owner:** Emmanuel Ogbu (Manny) — MSc Bioinformatics (Bradford), BSc Biomedical Science (MMU)
**Repo:** github.com/MannyMotion/Bioinformatics-Agent
**Local path:** `C:\Users\Invate\Downloads\Bioinformatics-Agent`
**Last updated:** 2026-07-07
**Latest committed version:** v0.5.1
**Working (uncommitted) version:** v0.6.0-dev — *the agentic reasoning layer*

---

## 1. What this project is (the one-paragraph pitch)

Upload any bioinformatics file. BioAgent auto-detects the type, **reasons about
what the data can answer**, proposes biological questions, runs the correct
analysis pipeline, generates interactive plots, and explains the results in
plain English using a local LLM (Ollama) grounded in a RAG knowledge base built
from Manny's MSc notes. The goal: let a biologist get from raw file to
interpreted result without knowing which tool to run or how to read the output.

**Why it exists:** Manny has the degrees (MSc Bioinformatics) but no industry
experience yet. This project is his portfolio proof that he can design agentic
systems, implement RAG, build full-stack apps, and reason about real biology.

---

## 2. Current architecture (as of the uncommitted working tree)

The system recently changed from a **fixed router** to a **two-step agentic flow**.
This is the most important recent change and is NOT yet reflected in the README.

```
Upload file
   │
   ▼
[detector.py]      → what file type is this? (FASTA/FASTQ/VCF/CSV, with confidence)
   │
   ▼
[analyst.py]       → PROFILE the data (dims, ID type, candidate groups, cleanliness)
   │                → PROPOSE biological questions the data can actually answer
   │                → INFER control vs treatment groups from column NAMES (not fixed prefixes)
   │                → PLAN which pipeline + params fit the chosen question
   │
   ▼  (user picks a proposed question OR types their own)
[pipelines/*]      → run the analysis (FASTA QC / RNA-seq DE / Variant annotation)
   │
   ▼
[explainer.py]     → RAG retrieve + Ollama(results + context) → plain-English interpretation
   │
   ▼
[job_store.py]     → persist the whole job to SQLite (jobs.db) so Q&A survives restarts
   │
   ▼
Frontend Q&A box   → /ask/{job_id} → follow-up questions answered from stored job
```

### API endpoints (`src/bioagent/api/main.py`)
- `GET  /health` — liveness check
- `POST /upload` — saves file, profiles it, returns proposed questions + inferred groups (fast, offline, no LLM)
- `POST /analyze/{job_id}` — runs the chosen pipeline for the picked question
- `POST /ask/{job_id}` — LLM Q&A against the stored job result (rate limit 20/min)

### Key modules
| File | Role | Status |
|------|------|--------|
| `agent/detector.py` | File-type detection w/ confidence | ✅ working, recently edited |
| `agent/analyst.py` | **NEW** — profile / propose / infer groups / plan | ⚠️ uncommitted |
| `agent/gene_ids.py` | **NEW** — Ensembl ID → HGNC symbol (curated map) | ⚠️ uncommitted |
| `agent/router.py` | Pipeline routing | ✅ working |
| `agent/explainer.py` | Ollama Q&A + RAG grounding (ANSI cleaned) | ✅ working |
| `pipelines/fasta_qc.py` | FASTA/FASTQ QC, Chart.js plots | ✅ working |
| `pipelines/rnaseq.py` | RNA-seq DE — now accepts explicit group cols | ✅ working, recently edited |
| `pipelines/variant_annotation.py` | VCF clinical annotation | ✅ working |
| `rag/*` | ChromaDB + sentence-transformers (all-MiniLM-L6-v2), 640 chunks | ✅ working |
| `api/job_store.py` | **NEW** — SQLite job persistence | ⚠️ uncommitted |
| `api/main.py` | FastAPI app, 4 endpoints, thread-pool for pipelines | ✅ working, heavily edited |
| `frontend/index.html` | Single-page UI: upload → questions → results → Q&A | ✅ working, heavily edited |
| `static/chart.umd.min.js` | **NEW** — Chart.js served locally (no CDN) | ⚠️ uncommitted |

---

## 3. Tech stack

- **Backend:** FastAPI + Uvicorn, Python 3.11
- **Frontend:** HTML/CSS/vanilla JS, Chart.js (served locally from `static/`)
- **LLM:** Ollama `llama3.2:3b` — 100% local, no API key
- **RAG:** ChromaDB + sentence-transformers (`all-MiniLM-L6-v2`), 640 chunks
- **Bioinformatics:** pandas, NumPy, SciPy, BioPython
- **Persistence:** SQLite (`jobs.db`)
- **Security:** slowapi rate limiting, 50 MB cap, extension whitelist, 500-char question cap
- **Env:** conda env `bioagent` → `C:\Users\Invate\anaconda3\envs\bioagent\python.exe`

---

## 4. How to run it

```bash
conda activate bioagent
cd C:\Users\Invate\Downloads\Bioinformatics-Agent
& "C:\Users\Invate\anaconda3\envs\bioagent\python.exe" -m uvicorn bioagent.api.main:app --host 0.0.0.0 --port 8000
# open http://localhost:8000/frontend/index.html
```
Tests: `& "C:\Users\Invate\anaconda3\envs\bioagent\python.exe" -m pytest tests/ -v`

**Sample data** (`data/sample/`): `test.fasta`, `test.vcf`, `counts.csv`,
`real_breast_cancer.csv`, `real_clinical_variants.vcf`, `ecoli_k12.fasta`.

---

## 5. Known issues / technical debt

- **PCA disabled** in `rnaseq.py` — `np.linalg.eigh` segfaults in the Windows
  uvicorn worker. Resolves on Linux deploy.
- **matplotlib removed** from all pipelines (same Windows segfault) — plots are
  now Chart.js HTML files written directly. The `*_plot_runner.py` files in
  `utils/` are UNUSED remnants of the subprocess approach.
- **DE method is a two-sample t-test on CPM.** Works, gives biologically
  sensible results on the test data, but is *not* what production RNA-seq uses
  (see §7 — this is the #1 credibility item for a bioinformatics interviewer).
- Ollama `3b` answers can be slightly repetitive — upgrade to 7b when hardware allows.
- ChromaDB telemetry warnings — cosmetic, ignore.

---

## 6. Validated results (real data, still true)

- **Breast cancer RNA-seq** — correctly surfaced ERBB2/HER2 (+4.05 log2FC),
  MKI67 (+4.84); FOXA1, BRCA1 down. Real clinical biomarkers.
- **Clinical VCF** — 8 pathogenic variants (BRCA1, BRCA2, TP53, KRAS, EGFR,
  PIK3CA, NRAS) flagged from ClinVar rsIDs.
- **E. coli K-12 FASTA** — mean GC 47.37% (published ~50.8%), no low-complexity.

---

## 7. Roadmap & priorities (ordered by job-market impact)

**Do next (highest leverage for landing a job):**
1. **Commit the agentic layer.** `analyst.py`, `gene_ids.py`, `job_store.py`,
   `static/`, and the main.py/frontend/rnaseq changes are unsaved. Tag v0.6.0.
2. **Upgrade the DE statistics** — swap the t-test for a proper method
   (pydeseq2 = DESeq2 in Python, or at minimum add multiple-testing correction
   with Benjamini-Hochberg FDR). This is the single biggest credibility fix.
3. **Deploy a live demo** (Render/Railway/HF Spaces) so applications link a URL,
   not a localhost screenshot. Re-enable PCA + matplotlib on the Linux box.
4. **Write the methods honestly** in the README — say what stats were used and
   why. Interviewers respect "t-test as a v1, DESeq2 planned" far more than
   silence.

**Later:**
- Proteomics pipeline (mass-spec CSV), metagenomics pipeline.
- User accounts + job history UI.
- FASTQ quality trimming.

---

## 8. Session log (newest first — append one block per session)

### 2026-07-07 (later) — Fable — Dashboard redesign
- Rebuilt `frontend/index.html` into a real product landing page while keeping
  ALL existing JS and element IDs intact (backend flow untouched). Verified
  every `getElementById` target still resolves.
- Added: sticky nav, hero with value prop + capability chips, upload card as the
  hero's primary action, a credibility strip (3 pipelines / 640 chunks / 100%
  local / ClinVar), a 4-step "how it works" flow, and three pipeline detail
  cards showing the actual biology. Genomic dark theme, dot-grid motif.
- Only JS logic change: `showAgent()`/`resetUI()` now toggle a `#landingView`
  wrapper (which contains the upload card) instead of `#uploadSection` alone.
- Rendered headless in Chrome to confirm layout. Looks clean and professional.
- **Next:** commit this + the agentic layer (still uncommitted), then the DE
  statistics upgrade (§7 item 2).

### 2026-07-07 — Fable
- Read the whole repo. Found the uncommitted **agentic layer** (`analyst.py` +
  friends) — the real differentiator, not yet committed or documented.
- Created this `PROJECT_MEMORY.md` as the new single source of truth.
- Gave Manny an honest assessment of the project + career direction (see chat):
  project is genuinely strong; #1 risk is being able to *explain* the code in
  interviews; #1 technical fix is upgrading the DE stats beyond a t-test.
- **Next session should start by:** deciding whether to (a) commit + tag the
  agentic layer, or (b) upgrade the DE statistics first.

### (historical) up to 2026-05-03 — Steve
- Built RAG, detector, router, 3 pipelines, FastAPI backend, frontend.
- Added Ollama Q&A (v0.5.0) and security hardening (v0.5.1).
- Validated all three pipelines on real biological datasets.
- Full detail preserved in `.claude/CONTINUITY.md`.

---

## 9. How to keep this memory current (the ritual)

At the **end of every working session**, update:
1. `Last updated` date + version fields at the top.
2. Any module status changes in §2.
3. New known issues in §5.
4. A new dated block at the top of §8 (what we did, what's next).

At the **start of every session**, read §1, §2, and the latest §8 block.
That's the whole "never forget where we are" system.
