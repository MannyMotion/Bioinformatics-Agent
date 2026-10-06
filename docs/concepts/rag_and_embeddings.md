# RAG and embeddings

**What it is:** Retrieval-Augmented Generation: look up relevant text first, then hand it to the LLM along with the question, so it answers using your material instead of guessing.

**Biology / intuition:** Imagine an open-book exam. The student (LLM) is better when handed the right pages (retrieved chunks from my MSc notes).

**How it works in BioAgent:**
1. **Chunking:** notes split into 500-character pieces with 100 characters of overlap so ideas aren't cut at boundaries ([ingestion.py](../../src/bioagent/rag/ingestion.py)).
2. **Embedding:** each chunk is turned into a vector of numbers by `all-MiniLM-L6-v2` (sentence-transformers) so that texts with similar meaning get nearby vectors ([embedder.py](../../src/bioagent/rag/embedder.py)).
3. **Store:** vectors go in ChromaDB on disk ([vector_store.py](../../src/bioagent/rag/vector_store.py)). About 640 chunks.
4. **Retrieve:** the question is embedded the same way, and the closest 2–3 chunks (by vector distance) are returned ([retriever.py](../../src/bioagent/rag/retriever.py)).
5. **Generate:** chunks plus results go into the Ollama prompt ([explainer.py](../../src/bioagent/agent/explainer.py)).

**When it's wrong / limits:**
- Retrieval can fetch irrelevant chunks, and the model can ignore or distort them. RAG reduces hallucination. It doesn't remove it. Our README's "grounded, not hallucinated" over-claims.
- The corpus is only as good as my notes.
- Our prompts include only summary statistics, not the actual gene list, so answers about specific genes come from the model's memory plus loosely related chunks.
- Short chunks can lose context.

**Likely interview questions**
1. *What is an embedding?* — A numeric vector representing a piece of text's meaning, so similar meanings sit close together and I can search by meaning, not keywords.
2. *Why RAG instead of fine-tuning?* — It's cheap, needs no training, and I can update the knowledge by adding documents. The cost is retrieval quality and prompt length.
