"""
embedder.py

Purpose: Convert text chunks into numerical vector representations (embeddings).
         Embeddings are how the RAG system understands meaning — two pieces of
         text with similar meaning will have similar vectors, even if they use
         different words.

         Example: "DNA sequence quality" and "nucleotide read accuracy" will
         have vectors close together in embedding space.

Inputs:  String or list of strings (text chunks)
Outputs: List of embedding vectors (list of floats)

Author:  Emmanuel Ogbu (Manny)
Date:    2026-04-23
"""

import logging
from sentence_transformers import SentenceTransformer

from bioagent.utils.logger import get_logger

logger = get_logger(__name__)

# Model choice: all-MiniLM-L6-v2 is a lightweight but powerful embedding model.
# 80MB on disk, runs on CPU, 384-dimensional output vectors.
# Benchmarked well on semantic similarity tasks.
# Source: https://huggingface.co/sentence-transformers/all-MiniLM-L6-v2
MODEL_NAME = "all-MiniLM-L6-v2"


class BioEmbedder:
    """
    Wraps sentence-transformers to produce embeddings for the RAG system.

    We use a class here (not just a function) because loading the model
    is expensive — ~2 seconds on first load. By storing it as an instance
    attribute, we load it once and reuse it for all subsequent calls.
    """

    def __init__(self) -> None:
        """Load the embedding model into memory on initialisation."""
        logger.info(f"Loading embedding model: {MODEL_NAME}")

        try:
            # SentenceTransformer downloads the model on first use (~80MB)
            # and caches it locally — subsequent loads are instant
            self.model = SentenceTransformer(MODEL_NAME)
            logger.info("Embedding model loaded successfully.")
        except Exception as e:
            logger.error(f"Failed to load embedding model: {e}")
            raise

    def embed(self, texts: list[str]) -> list[list[float]]:
        """
        Convert a list of text strings into embedding vectors.

        Args:
            texts: List of text strings to embed. Can be single sentences
                   or longer paragraphs — model handles both.

        Returns:
            List of embedding vectors. Each vector is a list of 384 floats.
            The order matches the input list exactly.

        Raises:
            ValueError: If texts list is empty.
        """
        if not texts:
            raise ValueError("Cannot embed an empty list of texts.")

        logger.info(f"Embedding {len(texts)} text chunk(s)...")

        # encode() returns a numpy array of shape (n_texts, 384)
        # convert_to_numpy=True is default but explicit here for clarity
        embeddings = self.model.encode(
            texts,
            convert_to_numpy=True,   # returns numpy array, not torch tensor
            show_progress_bar=False  # suppress progress bar in production
        )

        # Convert numpy array to plain Python list — ChromaDB expects lists
        return embeddings.tolist()

    def embed_single(self, text: str) -> list[float]:
        """
        Convenience method to embed a single string.

        Args:
            text: Single text string to embed.

        Returns:
            Single embedding vector as list of 384 floats.
        """
        return self.embed([text])[0]