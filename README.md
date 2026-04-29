# 🧬 BioAgent — Agentic Bioinformatics Analysis System

![Version](https://img.shields.io/badge/version-0.4.0-blue)
![Python](https://img.shields.io/badge/python-3.11-green)
![License](https://img.shields.io/badge/license-MIT-orange)
![Tests](https://img.shields.io/badge/tests-23%20passing-brightgreen)

> An AI-powered bioinformatics platform that automatically detects uploaded data files, selects the appropriate analysis pipeline, executes the full workflow, and explains every decision in plain English — backed by a RAG knowledge base trained on MSc-level bioinformatics notes.

---

## 🎯 What BioAgent Does

A researcher uploads any bioinformatics file. BioAgent:

1. **Auto-detects** the file type (FASTA, FASTQ, VCF, CSV) with confidence scoring
2. **Selects** the correct pipeline and explains why in plain English
3. **Executes** the full analysis automatically
4. **Generates** publication-quality visualisations
5. **Interprets** results using a RAG knowledge base of bioinformatics literature

No configuration. No pipeline selection. Just upload and analyse.

---

## 🏗️ Architecture

┌─────────────────────────────────────────────────────┐
│                   Web Frontend                       │
│         (Drag & drop upload, results display)        │
└─────────────────────┬───────────────────────────────┘
│ HTTP
┌─────────────────────▼───────────────────────────────┐
│                  FastAPI Backend                     │
│              POST /upload  GET /health               │
└──────────┬──────────────────────────────────────────┘
│
┌──────────▼──────────┐    ┌─────────────────────────┐
│   File Detector     │    │     RAG Knowledge Base   │
│  (auto-detects type)│    │  640 chunks from MSc     │
└──────────┬──────────┘    │  bioinformatics notes    │
│               │  ChromaDB + sentence-    │
┌──────────▼──────────┐    │  transformers            │
│   Pipeline Router   │◄───┘                          │
│  (selects pipeline) │                               │
└──────────┬──────────┘                               │
│                                          │
┌──────┴────────────────┐                        │
│                       │                        │
┌───▼────┐  ┌──────────┐  ┌▼──────────┐            │
│FASTA QC│  │ RNA-seq  │  │  Variant  │            │
│Pipeline│  │ Pipeline │  │Annotation │            │
└───┬────┘  └────┬─────┘  └─────┬─────┘            │
│            │              │                   │
└────────────┴──────────────┘                   │
│                                  │
┌────────────▼──────────────┐                  │
│   Visualisation Engine    │                  │
│ Matplotlib + Seaborn +    │                  │
│ Publication-quality plots │                  │
└───────────────────────────┘                  │

---

## 🔬 Pipelines

### Pipeline 1 — FASTA Quality Control
**Triggered by:** `.fasta`, `.fa`, `.fna`, `.fastq` files

Analyses:
- GC content per sequence (flags AT-rich or GC-rich outliers)
- Sequence length distribution
- Nucleotide composition (A/T/G/C/N ratios)
- Low complexity sequence detection (k-mer diversity)

Outputs:
- GC content distribution histogram
- Sequence length distribution plot
- Nucleotide composition bar chart
- QC report with biological interpretation

---

### Pipeline 2 — RNA-seq Differential Expression
**Triggered by:** `.csv`, `.tsv` gene count matrices

Analyses:
- CPM normalisation (removes library size bias)
- Differential expression (t-test, log2 fold change)
- Significance filtering (p < 0.05, |log2FC| ≥ 1.0)

Outputs:
- Volcano plot (fold change vs significance, genes labelled)
- PCA plot (sample clustering by condition)
- Heatmap (top differentially expressed genes)
- Full DE results with biological interpretation

---

### Pipeline 3 — Variant Annotation
**Triggered by:** `.vcf` files

Analyses:
- VCF parsing (CHROM, POS, REF, ALT, QUAL, INFO fields)
- Quality filtering (PHRED score threshold)
- Clinical annotation (gene name, consequence, condition)
- Pathogenicity classification

Outputs:
- Variant quality score chart
- Variants per gene bar chart
- Clinical significance pie chart
- Pathogenic variant report with gene-disease associations

---

## 🧠 RAG Knowledge Base

The system includes a semantic knowledge base built from MSc Bioinformatics lecture notes covering:

- Genomics and sequence analysis
- Transcriptomics and RNA-seq methods
- Proteomics
- Variant calling and annotation
- FastQC quality control
- Sequence alignment and phylogenetics
- BioPython workflows

**Technology:** ChromaDB vector database + sentence-transformers (`all-MiniLM-L6-v2`)
**Size:** 640 chunks, queryable in <2 seconds

---

## 🚀 Installation & Setup

### Prerequisites
- Python 3.11+
- Anaconda or Miniconda
- Git

### 1. Clone the repository
```bash
git clone https://github.com/MannyMotion/Bioinformatics-Agent.git
cd Bioinformatics-Agent
```

### 2. Create the conda environment
```bash
conda env create -f environment.yml
conda activate bioagent
```

### 3. Install the package
```bash
pip install -e .
```

### 4. Start the backend server
```bash
python -m uvicorn bioagent.api.main:app --reload --host 0.0.0.0 --port 8000
```

### 5. Open the frontend
Open `frontend/index.html` in your browser. You should see **"Server Online"** in the top right.

---

## 📁 Project Structure

bioinformatics-agent/
├── src/bioagent/
│   ├── agent/
│   │   ├── detector.py          # File type auto-detection
│   │   └── router.py            # Pipeline routing logic
│   ├── api/
│   │   └── main.py              # FastAPI REST backend
│   ├── parsers/
│   │   └── fasta_parser.py      # FASTA file parser
│   ├── pipelines/
│   │   ├── fasta_qc.py          # FASTA QC pipeline
│   │   ├── rnaseq.py            # RNA-seq pipeline
│   │   └── variant_annotation.py # Variant annotation pipeline
│   ├── rag/
│   │   ├── embedder.py          # Text embedding (sentence-transformers)
│   │   ├── ingestion.py         # Document chunking and indexing
│   │   ├── retriever.py         # Semantic similarity search
│   │   └── vector_store.py      # ChromaDB interface
│   └── utils/
│       └── logger.py            # Centralised logging
├── frontend/
│   └── index.html               # Web interface
├── data/
│   ├── sample/                  # Test files (FASTA, VCF, CSV)
│   └── knowledge/               # Ingested lecture notes
├── tests/                       # 23 unit tests
├── environment.yml              # Conda environment
└── requirements.txt             # Python dependencies

---

## 🧪 Running Tests

```bash
python -m pytest tests/ -v
```

Expected output: **23 passed**

---

## 🛠️ Tech Stack

| Layer | Technology |
|-------|-----------|
| Frontend | HTML5, CSS3, Vanilla JS |
| Backend | FastAPI, Uvicorn |
| RAG | ChromaDB, sentence-transformers |
| Bioinformatics | BioPython, pandas, scipy |
| Visualisation | Matplotlib, Seaborn |
| ML | scikit-learn (PCA) |
| Language | Python 3.11 |
| Version Control | Git, GitHub |

---

## 📊 Test Results

tests/test_detector.py          4 passed
tests/test_fasta_qc.py          5 passed
tests/test_retriever.py         4 passed
tests/test_rnaseq.py            5 passed
tests/test_variant_annotation.py 5 passed
─────────────────────────────────────────
TOTAL                          23 passed

---

## 🗺️ Roadmap

### v0.5.0 — LLM Integration (Ollama)
- Integrate Ollama (llama3.2:3b) as reasoning engine
- Natural language explanations of every analysis step
- Interactive Q&A about results
- Pipeline decision explanations grounded in RAG context

### v0.6.0 — Additional Pipelines
- Metagenomics pipeline (16S rRNA analysis)
- Protein structure prediction pipeline
- Epigenomics (ATAC-seq, ChIP-seq)

### v1.0.0 — SaaS Launch
- User authentication
- Subscription tiers (Free / Pro / Enterprise)
- Cloud deployment (AWS/GCP)
- API access for programmatic use
- Team collaboration features

---

## 👤 Author

**Emmanuel Ogbu (Manny)**
MSc Bioinformatics — University of Bradford
BSc Biomedical Science — Manchester Metropolitan University

GitHub: [@MannyMotion](https://github.com/MannyMotion)

---

## 📚 References

- Andrews S. (2010). FastQC: A Quality Control Tool for High Throughput Sequence Data. Babraham Bioinformatics.
- Love MI, Huber W, Anders S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15:550. https://doi.org/10.1186/s13059-014-0550-8
- McLaren W, et al. (2016). The Ensembl Variant Effect Predictor. *Genome Biology*, 17:122. https://doi.org/10.1186/s13059-016-0974-4
- Kim D, et al. (2019). Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. *Nature Biotechnology*, 37:907–915. https://doi.org/10.1038/s41587-019-0201-4