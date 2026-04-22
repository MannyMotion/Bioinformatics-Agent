"""
ingestion.py

Purpose: Ingest documents into the RAG knowledge base.
         Takes raw text files or PDFs, splits them into chunks,
         embeds each chunk, and stores them in ChromaDB.

         This is how the system "learns" from your lecture notes,
         research papers, and bioinformatics protocols.

Inputs:  Path to a text file or PDF
Outputs: Chunks stored in ChromaDB, ready for retrieval

Author:  Emmanuel Ogbu (Manny)
Date:    2026-04-23
"""

import re
import hashlib
from pathlib import Path

from bioagent.rag.embedder import BioEmbedder
from bioagent.rag.vector_store import BioVectorStore
from bioagent.utils.logger import get_logger

logger = get_logger(__name__)

# Chunk size in characters — how much text per chunk.
# Too large: chunks are too general, retrieval is imprecise.
# Too small: chunks lose context, retrieval misses meaning.
# 500 characters (~80-100 words) is a good starting point.
CHUNK_SIZE = 500

# Overlap between chunks — ensures context isn't lost at boundaries.
# If chunk 1 ends mid-sentence, chunk 2 starts 100 chars back to capture it.
CHUNK_OVERLAP = 100


def ingest_text_file(file_path: str | Path) -> int:
    """
    Ingest a plain text file into the RAG knowledge base.

    Reads the file, splits into overlapping chunks, embeds each chunk,
    and stores in ChromaDB. Returns the number of chunks stored.

    Args:
        file_path: Path to the .txt file to ingest.

    Returns:
        Number of chunks successfully stored.

    Raises:
        FileNotFoundError: If the file does not exist.
    """
    path = Path(file_path)

    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")

    logger.info(f"Ingesting file: {path.name}")

    # Read the full text
    try:
        text = path.read_text(encoding="utf-8")
    except Exception as e:
        logger.error(f"Failed to read file {path}: {e}")
        raise

    # Split into chunks
    chunks = _chunk_text(text, CHUNK_SIZE, CHUNK_OVERLAP)
    logger.info(f"Split into {len(chunks)} chunks.")

    if not chunks:
        logger.warning(f"No chunks generated from {path.name} — file may be empty.")
        return 0

    # Generate unique IDs for each chunk using file name + chunk index
    # We use a hash of the content to avoid duplicate storage
    ids = [_generate_chunk_id(path.name, i, chunk) for i, chunk in enumerate(chunks)]

    # Attach metadata so we know where each chunk came from
    metadatas = [
        {
            "source": path.name,
            "chunk_index": i,
            "total_chunks": len(chunks)
        }
        for i in range(len(chunks))
    ]

    # Embed all chunks in one batch — more efficient than one at a time
    embedder = BioEmbedder()
    embeddings = embedder.embed(chunks)

    # Store in ChromaDB
    store = BioVectorStore()

    # Filter out chunks that are already stored (avoid duplicates)
    new_chunks, new_embeddings, new_ids, new_metadatas = _filter_existing(
        store, chunks, embeddings, ids, metadatas
    )

    if not new_chunks:
        logger.info(f"All chunks from {path.name} already in knowledge base.")
        return 0

    store.add_documents(
        documents=new_chunks,
        embeddings=new_embeddings,
        ids=new_ids,
        metadatas=new_metadatas
    )

    logger.info(f"Ingested {len(new_chunks)} new chunks from {path.name}.")
    return len(new_chunks)


def _chunk_text(text: str, chunk_size: int, overlap: int) -> list[str]:
    """
    Split text into overlapping chunks of approximately chunk_size characters.

    Overlap ensures that sentences split across chunk boundaries are still
    captured in at least one complete chunk.

    Args:
        text: Full document text.
        chunk_size: Target size of each chunk in characters.
        overlap: Number of characters to overlap between chunks.

    Returns:
        List of text chunks.
    """
    # Clean up excessive whitespace — normalise the text first
    text = re.sub(r"\s+", " ", text).strip()

    if len(text) <= chunk_size:
        # Document is small enough to be a single chunk
        return [text]

    chunks = []
    start = 0

    while start < len(text):
        end = start + chunk_size

        # Try to break at a sentence boundary (". ") rather than mid-word
        if end < len(text):
            # Look for the last sentence end within the chunk
            boundary = text.rfind(". ", start, end)
            if boundary != -1 and boundary > start:
                end = boundary + 1  # include the period

        chunk = text[start:end].strip()
        if chunk:
            chunks.append(chunk)

        # Move forward by chunk_size minus overlap
        start += chunk_size - overlap

    return chunks


def _generate_chunk_id(filename: str, index: int, content: str) -> str:
    """
    Generate a unique, deterministic ID for a chunk.

    Uses MD5 hash of filename + index + first 50 chars of content.
    Deterministic means the same chunk always gets the same ID —
    this is how we detect and skip duplicates.

    Args:
        filename: Source file name.
        index: Chunk index within the file.
        content: The chunk text.

    Returns:
        Unique string ID.
    """
    raw = f"{filename}_{index}_{content[:50]}"
    # MD5 is fine here — we just need uniqueness, not cryptographic security
    return hashlib.md5(raw.encode()).hexdigest()


def _filter_existing(
    store: BioVectorStore,
    chunks: list[str],
    embeddings: list[list[float]],
    ids: list[str],
    metadatas: list[dict]
) -> tuple[list, list, list, list]:
    """
    Filter out chunks that are already stored in ChromaDB.

    Args:
        store: The BioVectorStore instance.
        chunks, embeddings, ids, metadatas: Parallel lists to filter.

    Returns:
        Filtered versions of all four lists containing only new chunks.
    """
    # Get all existing IDs from ChromaDB
    try:
        existing = store.collection.get(ids=ids)
        existing_ids = set(existing["ids"])
    except Exception:
        # If get() fails (e.g. empty collection), assume nothing exists
        existing_ids = set()

    # Keep only chunks whose IDs aren't already stored
    new_chunks, new_embeddings, new_ids, new_metadatas = [], [], [], []
    for chunk, emb, id_, meta in zip(chunks, embeddings, ids, metadatas):
        if id_ not in existing_ids:
            new_chunks.append(chunk)
            new_embeddings.append(emb)
            new_ids.append(id_)
            new_metadatas.append(meta)

    return new_chunks, new_embeddings, new_ids, new_metadatas