# BioAgent — Project Continuity Document

**Last Updated:** April 28, 2026
**Current Version:** v0.4.0
**Project Status:** Demo-ready, all three pipelines functional

---

## 🎯 Project Summary

BioAgent is an agentic bioinformatics analysis system that:
1. Auto-detects uploaded bioinformatics files (FASTA, FASTQ, VCF, CSV)
2. Routes to the correct pipeline
3. Executes full analysis with publication-quality plots
4. Explains results using a RAG knowledge base (640 chunks from MSc notes)
5. Displays everything in a web browser

**Current milestone:** Full-stack demo complete (frontend + backend + all 3 pipelines).
**Next milestone:** Ollama LLM integration (v0.5.0).

---

## 📦 What's Built (v0.4.0)

### Layer 1 — RAG Knowledge System
- **ChromaDB** vector database (local, free)
- **sentence-transformers** (`all-MiniLM-L6-v2`) for embeddings
- **640 chunks** ingested from MSc Bioinformatics notes (Sem1, Sem2, Sem3)
- Queryable in <2 seconds

**Key files:**
- `src/bioagent/rag/embedder.py` — text to vectors
- `src/bioagent/rag/vector_store.py` — ChromaDB interface
- `src/bioagent/rag/ingestion.py` — document chunking
- `src/bioagent/rag/retriever.py` — semantic search

---

### Layer 2 — File Detection & Routing
- **Auto-detector** (`src/bioagent/agent/detector.py`)
  - FASTA/FASTQ detection (95%+ confidence)
  - VCF detection (99%+ confidence)
  - CSV/TSV detection (80%+ confidence)
  - Returns `DetectionResult` with confidence and explanation

- **Router** (`src/bioagent/agent/router.py`)
  - Routes detected files to correct pipeline
  - Provides human-readable reasoning for each routing decision

---

### Layer 3 — Three Bioinformatics Pipelines

#### Pipeline 1: FASTA Quality Control (`src/bioagent/pipelines/fasta_qc.py`)
- **Input:** FASTA or FASTQ files
- **Analyses:**
  - GC content per sequence
  - Sequence length distribution
  - Nucleotide composition
  - Low complexity detection (k-mer diversity)
- **Outputs:** 3 plots (GC histogram, length distribution, composition bar chart) + interpretation

#### Pipeline 2: RNA-seq Differential Expression (`src/bioagent/pipelines/rnaseq.py`)
- **Input:** Gene count matrix CSV (rows=genes, cols=samples)
- **Analyses:**
  - CPM normalisation
  - T-test differential expression
  - Significance filtering (p<0.05, |log2FC|≥1.0)
- **Outputs:** Volcano plot, PCA plot, heatmap + full DE report

#### Pipeline 3: Variant Annotation (`src/bioagent/pipelines/variant_annotation.py`)
- **Input:** VCF files
- **Analyses:**
  - VCF parsing
  - Quality filtering (PHRED threshold)
  - Clinical annotation (gene, consequence, pathogenicity)
- **Outputs:** Quality chart, genes-per-variant chart, significance pie chart + pathogenic report

---

### Layer 4 — FastAPI Backend & Web Frontend

**Backend** (`src/bioagent/api/main.py`):
- `POST /upload` — accepts file + optional RNA-seq labels, returns analysis results
- `GET /health` — server status
- `GET /analyse/{job_id}` — retrieve previous results
- `GET /jobs` — list all completed jobs
- CORS enabled for frontend
- `/outputs` mounted as static file server for plot display

**Frontend** (`frontend/index.html`):
- Drag-and-drop file upload
- Auto-shows RNA-seq options for CSV files
- Animated progress steps during analysis
- Results display with stats cards, reasoning box, warnings, 3 plots, full interpretation
- "Run Another Analysis" button to reset

---

## 🗂️ File Structure