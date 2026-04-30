# 🛠️ Skills — Reusable Workflows

## Skill 1: How to Build a New Bioinformatics Pipeline

**Trigger:** "Let's build the metagenomics pipeline" or similar

**Process:**
1. Start with biology context — what does this pipeline analyse? Why?
2. Cite source: PubMed paper, official docs, Nature Methods
3. Define inputs (file format) and outputs (visualisations + stats)
4. Create Python module: `src/bioagent/pipelines/{name}.py`
5. Structure:
   - Top docstring: purpose, inputs, outputs, author, date
   - @dataclass for result object
   - Main function that orchestrates the pipeline
   - Private helper functions for each major step
6. Generate 2-3 publication-quality plots using matplotlib/seaborn
7. Write tests in `tests/test_{name}.py` (minimum 5 tests)
8. Integrate with router (`src/bioagent/agent/router.py`)
9. Test end-to-end: upload file → detect → route → analyse → display in browser
10. Commit with: `feat: add {name} pipeline`
11. Tag version (v0.x.0) when pipeline is complete

## Skill 2: How to Add New Knowledge to RAG

**Trigger:** "Let's ingest the metagenomics papers" or similar

**Process:**
1. Collect source documents (PDFs, lecture notes, papers)
2. Extract text: copy from PDF or Google Docs → save as `.txt` in `data/knowledge/`
3. Run ingestion:
```bash
   python -c "from bioagent.rag.ingestion import ingest_text_file; ingest_text_file('data/knowledge/myfile.txt')"
```
4. Verify chunks loaded: check ChromaDB count increased
5. Test retrieval: ask a question related to the new knowledge
6. Commit with: `chore: ingest {name} into RAG knowledge base`

## Skill 3: How to Deploy Locally + Test End-to-End

**Trigger:** "Let's test the full system" or new feature added

**Process:**
1. Activate conda: `conda activate bioagent`
2. Start backend in terminal 1: `python -m uvicorn bioagent.api.main:app --reload --host 0.0.0.0 --port 8000`
3. Wait for "Application startup complete"
4. Open frontend in browser: `frontend/index.html`
5. Should see "Server Online" green dot
6. Test each pipeline:
   - FASTA: upload `data/sample/test.fasta` → should show 2 sequences, 45.84% GC
   - VCF: upload `data/sample/test.vcf` → should show 8 pathogenic variants
   - RNA-seq: upload `data/sample/counts.csv` with control=healthy, treatment=cancer → should show 7 up, 3 down
7. Verify plots render in browser
8. Check server terminal for no errors
9. Run tests: `pytest tests/ -v`
10. If all pass, proceed to commit

## Skill 4: How to Write Tests (Minimum Viable)

**Trigger:** New module added

**Process:**
1. Create `tests/test_{module}.py`
2. Import the function/class you're testing
3. Create a fixture if needed (sample data file path)
4. Write minimum 5 tests:
   - Test happy path (correct input → expected output)
   - Test error handling (missing file → raises FileNotFoundError)
   - Test edge case (empty file, malformed data)
   - Test output type (returns correct dataclass, not dict)
   - Test integration (if applicable)
5. Each test: `def test_descriptive_name(fixture=None) -> None:`
6. Use `assert` statements, not print statements
7. Run: `pytest tests/test_{module}.py -v`
8. All tests must pass before commit

## Skill 5: How to Make a Clean Commit

**Trigger:** Feature complete

**Process:**
1. Run tests: `pytest tests/ -v`
2. Verify all pass
3. In GitHub Desktop, review all changes in "Changes" tab
4. Uncheck files that shouldn't be committed (data/*, outputs/*, uploads/*)
5. Write summary (under 72 chars, imperative mood):
   - `feat: add metagenomics pipeline`
   - `fix: handle empty FASTQ files`
   - `docs: update README with new section`
6. Write description (WHY, not what):