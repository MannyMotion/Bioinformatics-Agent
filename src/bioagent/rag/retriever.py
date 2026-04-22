"""
retriever.py

Purpose: Query interface for the RAG knowledge base.
         This is the module the bioinformatics analyzer calls when it
         needs context. It takes a plain English question, converts it
         to a vector, searches ChromaDB, and returns the most relevant
         knowledge chunks with their sources.

         This completes the RAG layer:
         ingestion.py → stores knowledge
         retriever.py → retrieves knowledge

Inputs:  A plain English question (string)
Outputs: List of relevant text chunks with metadata and relevance scores

Author:  Emmanuel Ogbu (Manny)
Date:    2026-04-23
"""

from dataclasses import dataclass
from bioagent.rag.embedder import BioEmbedder
from bioagent.rag.vector_store import BioVectorStore
from bioagent.utils.logger import get_logger

logger = get_logger(__name__)

# Maximum distance threshold — results further than this are too irrelevant
# to be useful. Tune this as the knowledge base grows.
MAX_DISTANCE = 1.2


@dataclass
class RetrievalResult:
    """
    A single retrieved knowledge chunk with its context.

    Using a dataclass so downstream modules get a clean,
    type-hinted object instead of raw dictionaries.
    """
    text: str           # the actual knowledge chunk
    source: str         # which file it came from
    chunk_index: int    # position within that file
    distance: float     # similarity score (lower = more relevant)
    relevance: float    # human-readable score (higher = more relevant)


class BioRetriever:
    """
    Retrieves relevant knowledge chunks from the RAG knowledge base.

    Initialises the embedder and vector store once on startup,
    then handles any number of queries efficiently.
    """

    def __init__(self) -> None:
        """Initialise embedder and vector store."""
        logger.info("Initialising BioRetriever...")
        # Load once — reused for every query
        self.embedder = BioEmbedder()
        self.store = BioVectorStore()
        logger.info(
            f"Retriever ready. Knowledge base contains "
            f"{self.store.count()} chunks."
        )

    def retrieve(
        self,
        question: str,
        n_results: int = 3
    ) -> list[RetrievalResult]:
        """
        Retrieve the most relevant knowledge chunks for a question.

        Args:
            question: Plain English question or bioinformatics query.
            n_results: Number of results to return (default 3).

        Returns:
            List of RetrievalResult objects, sorted by relevance.
            Empty list if no relevant results found above threshold.

        Raises:
            ValueError: If question is empty.
        """
        if not question.strip():
            raise ValueError("Question cannot be empty.")

        logger.info(f"Retrieving context for: '{question}'")

        # Convert question to vector
        query_vector = self.embedder.embed_single(question)

        # Search ChromaDB for nearest neighbours
        raw_results = self.store.query(query_vector, n_results=n_results)

        # Filter and package results
        results = []
        for r in raw_results:
            distance = r["distance"]

            # Skip results that are too dissimilar to be useful
            if distance > MAX_DISTANCE:
                logger.debug(
                    f"Skipping chunk (distance {distance:.4f} > {MAX_DISTANCE})"
                )
                continue

            # Convert distance to a 0-100 relevance score for readability
            # distance=0 → relevance=100, distance=MAX_DISTANCE → relevance=0
            relevance = round((1 - distance / MAX_DISTANCE) * 100, 1)

            results.append(RetrievalResult(
                text=r["text"],
                source=r["metadata"].get("source", "unknown"),
                chunk_index=r["metadata"].get("chunk_index", -1),
                distance=distance,
                relevance=relevance
            ))

        logger.info(f"Retrieved {len(results)} relevant chunk(s).")
        return results

    def retrieve_as_context(self, question: str, n_results: int = 3) -> str:
        """
        Retrieve relevant chunks and format them as a single context string.

        This is what gets passed to the LLM (Ollama) when we need the AI
        to answer a question grounded in our knowledge base.

        Args:
            question: The query to retrieve context for.
            n_results: Number of chunks to include in context.

        Returns:
            Formatted string with all relevant chunks and their sources.
        """
        results = self.retrieve(question, n_results)

        if not results:
            return "No relevant knowledge found in the knowledge base."

        # Format chunks into a clean context block
        context_parts = []
        for i, r in enumerate(results):
            context_parts.append(
                f"[Source {i+1}: {r.source} | Relevance: {r.relevance}%]\n"
                f"{r.text}"
            )

        return "\n\n".join(context_parts)