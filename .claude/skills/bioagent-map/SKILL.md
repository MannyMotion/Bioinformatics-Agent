---
name: bioagent-map
description: Codebase map for BioAgent. Use before reading or changing any BioAgent code, when asked "where is X", "how does a request flow", "how do I run or test this", or when a lesson needs to locate a module.
---

# BioAgent codebase map

Last verified: 2026-10-06 against git commit 7099cdf. Re-verify if files have moved. For rules and environment details, see [CLAUDE.md](../../../CLAUDE.md). For current state and known bugs, see [PROJECT_MEMORY.md](../../../PROJECT_MEMORY.md).

## Request flow
```
browser (frontend/index.html, vanilla JS, API = http://localhost:8000)
  POST /upload      -> api/main.py upload_file
        detector.detect_file_type   (first 20 lines -> FASTA/FASTQ/VCF/GFF/CSV/TSV)
        analyst.profile_file        (dims, ID type, infer groups, propose questions)
        job_store.save_job          (SQLite jobs.db, status "profiled")
  POST /analyze/{id} -> analyze_job
        analyst.plan_analysis       (pick pipeline + params)
        router.route_file           (re-detects, runs a pipeline)
          pipelines/fasta_qc.py | rnaseq.py | variant_annotation.py
        explainer.explain_results   (RAG retrieve + Ollama CLI; template fallback)
        job_store.save_job          (status "complete")
  POST /ask/{id}     -> ask_question: explainer.answer_question (RAG + Ollama)
  GET  /health
Plots: standalone Chart.js HTML written to outputs/, served at /outputs and /static/chart.umd.min.js.
```

## Modules (all under `src/bioagent/`)
| File | Job |
|---|---|
| `agent/detector.py` | Content-based file-type detection with confidence |
| `agent/analyst.py` | Profile data, infer control/treatment groups (keyword vocab), propose questions, plan. Rule-based, no LLM |
| `agent/gene_ids.py` | Ensembl ID detection and a small Ensembl-to-symbol map |
| `agent/router.py` | Route a file to a pipeline; returns (result, RoutingDecision) |
| `agent/explainer.py` | Build prompts, RAG context, call `ollama run llama3.2:3b` by subprocess |
| `pipelines/fasta_qc.py` | GC, lengths, composition, low-complexity (4-mer ratio), warnings, plots |
| `pipelines/rnaseq.py` | Filter, CPM, per-gene t-test + log2FC, volcano and heatmap, warnings |
| `pipelines/variant_annotation.py` | VCF parse, QUAL/FILTER filter, `KNOWN_VARIANTS` dict lookup, plots |
| `parsers/fasta_parser.py` | `parse_fasta` returning `FastaRecord`s (own parser, not Biopython) |
| `rag/` | `ingestion` (chunk 500/100), `embedder` (MiniLM), `vector_store` (ChromaDB at repo-root `chroma_db/`), `retriever` |
| `api/main.py` | FastAPI app, upload security checks, thread pool, rate limits (slowapi) |
| `api/job_store.py` | SQLite table `jobs(job_id, data JSON text, created_at)` at `jobs.db` |
| `utils/logger.py` | Shared logger. `utils/*_plot_runner.py` are UNUSED leftovers |

## Data and outputs
- Sample inputs: `data/sample/` (small, toy-sized; see [DATA_SOURCES.md](../../../DATA_SOURCES.md) for provenance). RAG source notes: `data/knowledge/`. Unversioned big data goes in `data/raw/` (ignored).
- Runtime (git-ignored): `uploads/`, `outputs/`, `jobs.db`, `chroma_db/`, `bioagent.log`.
- Tests: `tests/`. All but `test_api.py` run offline (see the live-server caveat in CLAUDE.md).

## Run and test
Commands live in the Environment section of [CLAUDE.md](../../../CLAUDE.md). Use the conda env python, not the system one. Run from the repo root.

## Gotchas
- RNA-seq: `_differential_expression` is given raw `filtered_df`, not CPM. PCA is disabled (Windows segfault in eigh). matplotlib imports remain at the top of `fasta_qc.py`/`variant_annotation.py`.
- Low-complexity flag fires for any sequence over about 860 bp (see [docs/concepts/low_complexity_sequences.md](../../../docs/concepts/low_complexity_sequences.md)).
- `/ask` passes only summary stats and the interpretation text to the LLM, not per-gene results.
- `main.py` version strings say 0.5.0. README does not yet describe the analyst layer.
- Job storage: `job_store.save_job` is INSERT OR REPLACE, so `/analyze` overwrites the `/upload` record (profile is copied across).
- `tests/` write plots into `outputs/` as a side effect.
- Some `outputs/*.png` files are tracked in git despite `outputs/` being ignored.
