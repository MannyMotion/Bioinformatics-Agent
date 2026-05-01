# Skill: Building Bioinformatics Pipelines

## How Manny Builds Pipelines (The Process)

### Step 1: Define the Biology First
Before any code, answer:
- What biological question does this data answer?
- What are the industry-standard tools? (cite sources)
- What are the expected inputs and outputs?

Example: For RNA-seq, the question is "which genes change expression between conditions?" Industry standard is DESeq2 (Love et al., 2014). Inputs: count matrix. Outputs: log2FC, p-values, volcano plot.

### Step 2: Design the Pipeline Function Signature
```python
def run_[pipeline_name](file_path, output_dir, use_rag=True) -> [ResultDataclass]:
    """
    Docstring with:
    - Purpose (one sentence)
    - What it does (numbered steps)
    - Args with types
    - Returns with type
    - Raises with conditions
    """
```

### Step 3: Build Step by Step
1. Load/parse the data (wrapped in try/except)
2. Validate/filter (log decisions)
3. Normalise (if needed, with comments explaining why)
4. Analyse (main computation)
5. Visualise (plots with publication quality)
6. Interpret (basic text first, RAG context later)

### Step 4: Test Before Shipping
- Write a test that uploads sample data
- Verify stats are correct
- Verify plots exist and render
- Verify interpretation is sensible

### Step 5: Document Everything
- Top-of-file docstring with author + date
- Every function has a docstring
- Inline comments on non-trivial lines
- Cite papers where biology is non-obvious

## What a Good Pipeline Has

### Correct Biology
- Uses industry-standard methods
- Cites authoritative sources
- Produces results that match published examples

### Clean Code
- Type hints on all parameters
- Docstrings on all functions
- Logging at key steps (not print)
- Try/except around file I/O and computation
- Guard clauses (fail early)

### Publication-Quality Plots
- 10x6 inch figures minimum
- Seaborn styling (`sns.set_theme()`)
- Axis labels with units
- Legends with sample counts
- DPI 150 minimum
- White space, not busy

### User-Facing Output
- A dataclass with results
- A warnings list (what to watch for)
- An interpretation string (plain English)
- Plot file paths
- Confidence/quality scores where applicable

## The Three Pipelines Manny Built

### Pipeline 1: FASTA QC
**Biology:** Assess sequence data quality before downstream analysis
**Key Methods:**
- GC content (organism signature, flags contamination)
- Length distribution (assembly coverage indicator)
- K-mer complexity (flags homopolymer runs)

**Source:** FastQC documentation (Andrews, Babraham Bioinformatics)

### Pipeline 2: RNA-seq DE
**Biology:** Find genes that change expression between conditions
**Key Methods:**
- CPM normalisation (removes library size bias)
- T-test (parametric test for means comparison)
- Volcano plot (2D significance view)

**Source:** Love et al. (2014) Genome Biology — DESeq2 paper

### Pipeline 3: Variant Annotation
**Biology:** Classify genetic variants by clinical impact
**Key Methods:**
- VCF parsing (standard genomics file format)
- Quality filtering (PHRED score > threshold)
- Gene annotation (which gene is affected)
- Pathogenicity lookup (known disease associations)

**Source:** McLaren et al. (2016) Genome Biology — VEP paper

## Next Pipeline Ideas (Not Yet Built)
- Metagenomics (16S rRNA taxonomy)
- Protein structure (AlphaFold integration)
- Epigenomics (ATAC-seq, ChIP-seq)