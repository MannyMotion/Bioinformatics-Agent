# 🧬 BioAgent — Agentic Bioinformatics Analysis System

> Upload any bioinformatics file. BioAgent auto-detects the data type, selects the correct analysis pipeline, runs the full workflow, and explains the results in plain English — powered by a local LLM and a RAG knowledge base.

![FASTA QC Pipeline](screenshots/fasta_qc.png)

---

## What It Does

Most bioinformatics tools require you to already know which pipeline to run, which parameters to set, and how to interpret the output. BioAgent removes that barrier.

You upload a file. The system does the rest:

1. **Auto-detects** the file type (FASTA, FASTQ, VCF, CSV, TSV)
2. **Selects** the correct pipeline and explains *why*
3. **Runs** the full analysis — QC, normalisation, differential expression, variant annotation
4. **Generates** interactive visualisations
5. **Interprets** the results in plain biological English
6. **Answers** follow-up questions via a local LLM (Ollama + RAG)

---

## Pipelines

### 🔬 FASTA/FASTQ — Sequence Quality Control
- GC content distribution
- Sequence length distribution
- Low complexity sequence detection
- Nucleotide composition analysis
- Biological interpretation of QC results

**Example result:** E. coli K-12 genome sequences — 47.37% GC content, median length 140bp, all QC metrics passing

![FASTA QC](screenshots/fasta_qc.png)

---

### 📊 RNA-seq — Differential Expression Analysis
- CPM normalisation
- Two-sample t-test for differential expression
- Volcano plot (log2 fold change vs significance)
- Heatmap of top differentially expressed genes
- Biological interpretation of gene findings

**Example result:** Breast cancer vs normal — identified ERBB2/HER2 (+4.05 log2FC), MKI67 (+4.84 log2FC) upregulated; FOXA1 (-2.55 log2FC), BRCA1 (-1.48 log2FC) downregulated. These are real clinically validated biomarkers.

![RNA-seq Analysis](screenshots/rnaseq.png)

---

### 🧬 Variant Annotation — Clinical VCF Analysis
- VCF parsing and quality filtering
- Clinical significance annotation (ClinVar rsIDs)
- Gene-level variant mapping
- Pathogenicity assessment
- Clinical alert system for pathogenic findings

**Example result:** 12 variants analysed — 8 pathogenic variants identified including BRCA1, BRCA2, TP53, KRAS, EGFR, PIK3CA

![Variant Annotation](screenshots/variant.png)

---

## AI-Powered Q&A

After every analysis, users can ask follow-up questions in plain English:

> *"What does ERBB2 upregulation mean clinically?"*
> *"Why is BRCA1 downregulated in cancer?"*
> *"What should I do next with these variant findings?"*

Ollama (llama3.2:3b) answers using the pipeline results + a RAG knowledge base built from bioinformatics literature — grounded, not hallucinated.

![Ollama Q&A](screenshots/qa.png)

---

## Architecture

```
┌─────────────────────────────────────────────────────────┐
│                    BioAgent System                       │
├─────────────────────────────────────────────────────────┤
│  Frontend (HTML/CSS/JS)                                  │
│  └── Drag-drop upload, interactive plots, Q&A box        │
├─────────────────────────────────────────────────────────┤
│  FastAPI Backend                                         │
│  ├── File Type Detector (FASTA/VCF/CSV — 95-99% conf.)  │
│  ├── Pipeline Router (selects correct analysis)          │
│  └── Security (rate limiting, file validation)           │
├─────────────────────────────────────────────────────────┤
│  Analysis Pipelines                                      │
│  ├── FASTA QC Pipeline                                   │
│  ├── RNA-seq Differential Expression Pipeline            │
│  └── Variant Annotation Pipeline                         │
├─────────────────────────────────────────────────────────┤
│  AI Layer                                                │
│  ├── RAG System (ChromaDB + sentence-transformers)       │
│  │   └── 640 chunks from bioinformatics literature       │
│  └── Ollama LLM (llama3.2:3b — runs 100% locally)       │
└─────────────────────────────────────────────────────────┘
```

---

## Tech Stack

| Layer | Technology |
|-------|-----------|
| Backend | FastAPI, Python 3.11 |
| Frontend | HTML, CSS, JavaScript, Chart.js |
| AI / LLM | Ollama (llama3.2:3b) — local, no API key |
| RAG | ChromaDB, sentence-transformers (all-MiniLM-L6-v2) |
| Data | pandas, NumPy, SciPy, BioPython |
| Security | slowapi (rate limiting), file validation |
| Version Control | Git / GitHub |

---

## Project Structure

```
Bioinformatics-Agent/
├── src/bioagent/
│   ├── agent/
│   │   ├── detector.py        # Auto file type detection
│   │   ├── router.py          # Pipeline selection logic
│   │   └── explainer.py       # Ollama Q&A integration
│   ├── pipelines/
│   │   ├── fasta_qc.py        # FASTA/FASTQ quality control
│   │   ├── rnaseq.py          # RNA-seq differential expression
│   │   └── variant_annotation.py  # VCF clinical annotation
│   ├── rag/
│   │   ├── embedder.py        # sentence-transformers embedding
│   │   ├── vector_store.py    # ChromaDB vector database
│   │   ├── ingestion.py       # Knowledge base ingestion
│   │   └── retriever.py       # RAG context retrieval
│   ├── api/
│   │   └── main.py            # FastAPI application
│   └── utils/
│       └── logger.py          # Structured logging
├── frontend/
│   └── index.html             # Single-page web interface
├── data/
│   ├── knowledge/             # RAG knowledge base documents
│   └── sample/                # Test datasets
├── tests/                     # 23 passing unit tests
├── scripts/                   # Utility scripts
└── requirements.txt
```

---

## Getting Started

### Prerequisites
- Python 3.11+
- Conda (recommended)
- [Ollama](https://ollama.com/download) installed

### Installation

```bash
# Clone the repository
git clone https://github.com/MannyMotion/Bioinformatics-Agent.git
cd Bioinformatics-Agent

# Create and activate conda environment
conda create -n bioagent python=3.11
conda activate bioagent

# Install dependencies
pip install -e .

# Pull the Ollama model (2GB download)
ollama pull llama3.2:3b
```

### Run the system

```bash
# Start the backend server
python -m uvicorn bioagent.api.main:app --host 0.0.0.0 --port 8000

# Open in browser
# http://localhost:8000/frontend/index.html
```

### Run tests

```bash
pytest tests/ -v
# Expected: 23 passed
```

---

## Supported File Formats

| Format | Pipeline | Use Case |
|--------|----------|----------|
| `.fasta`, `.fa`, `.fna` | FASTA QC | Genome sequences, assembled contigs |
| `.fastq` | FASTA QC | Raw sequencing reads |
| `.vcf` | Variant Annotation | SNPs, indels, clinical variants |
| `.csv`, `.tsv` | RNA-seq DE | Gene expression count matrices |

---

## Validated Results

System tested against real biological datasets:

- **Breast cancer RNA-seq** — correctly identified HER2/ERBB2 overexpression (+4.05 log2FC), a known therapeutic target in breast cancer
- **Clinical VCF** — detected 8 pathogenic variants including BRCA1, TP53, KRAS from ClinVar rsIDs
- **E. coli K-12 FASTA** — GC content 47.37% (published value: 50.8%), no low-complexity sequences detected

---

## Security

- Rate limiting: 10 uploads/minute per IP (slowapi)
- File size cap: 50MB maximum
- Extension whitelist: only known bioinformatics formats accepted
- Question length cap: 500 characters (prevents prompt injection)
- No API keys required — all AI runs locally via Ollama

---

## Version History

| Version | Description |
|---------|-------------|
| v0.1.0 | FASTA QC pipeline |
| v0.2.0 | RNA-seq differential expression pipeline |
| v0.3.0 | Variant annotation pipeline |
| v0.4.0 | Full-stack web interface |
| v0.5.0 | Ollama LLM integration + RAG Q&A |
| v0.5.1 | Security hardening (rate limiting, validation) |

---

## About

Built by **Emmanuel Ogbu (Manny)** — MSc Bioinformatics (University of Bradford), BSc Biomedical Science (Manchester Metropolitan University).

This project was built to demonstrate end-to-end agentic bioinformatics — combining classical computational biology pipelines with modern AI (RAG + LLM) to make genomic analysis accessible without requiring deep tool expertise.

**GitHub:** [github.com/MannyMotion](https://github.com/MannyMotion)

---

## Roadmap

- [ ] Deploy to cloud (Render/Railway) for live demo URL
- [ ] Re-enable PCA on Linux deployment
- [ ] FASTQ quality trimming (Trimmomatic integration)
- [ ] Metagenomics pipeline
- [ ] Proteomics pipeline (CSV mass spec data)
- [ ] User accounts and job history (SQLite)
- [ ] SaaS deployment

---

## Acknowledgements

- [FastAPI](https://fastapi.tiangolo.com/) — backend framework
- [Ollama](https://ollama.com/) — local LLM inference
- [ChromaDB](https://www.trychroma.com/) — vector database
- [Chart.js](https://www.chartjs.org/) — interactive visualisations
- [sentence-transformers](https://www.sbert.net/) — text embeddings
- BioPython, pandas, SciPy — bioinformatics data processing