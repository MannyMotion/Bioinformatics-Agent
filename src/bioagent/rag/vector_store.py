"""
vector_store.py

Purpose: Store and retrieve text embeddings using ChromaDB.
         ChromaDB is a local vector database — it stores text chunks
         alongside their embeddings and lets us find the most similar
         chunks to any query in milliseconds.

         Think of it as a search engine that understands meaning,
         not just keywords.

Inputs:  Text chunks + their embeddings (for storage)
         Query embedding (for retrieval)
Outputs: Most relevant text chunks for a given query

Author:  Emmanuel Ogbu (Manny)
Date:    2026-04-23
"""

from pathlib import Path
import chromadb
from chromadb.config import Settings

from bioagent.utils.logger import get_logger

logger = get_logger(__name__)

# ChromaDB will persist its database to this folder on disk
# This means knowledge survives between sessions — we don't re-embed every time
CHROMA_DB_PATH = Path(__file__).parent.parent.parent.parent / "chroma_db"

# Collection name — like a table in a regular database
COLLECTION_NAME = "bioagent_knowledge"


class BioVectorStore:
    """
    Manages storage and retrieval of text embeddings in ChromaDB.

    ChromaDB organises vectors into "collections" — similar to tables
    in SQL. We use one collection for all bioinformatics knowledge:
    lecture notes, papers, protocols.
    """

    def __init__(self) -> None:
        """Initialise ChromaDB client and get or create the collection."""
        logger.info(f"Initialising ChromaDB at: {CHROMA_DB_PATH}")

        try:
            # PersistentClient saves data to disk automatically
            # so knowledge persists between sessions
            self.client = chromadb.PersistentClient(
                path=str(CHROMA_DB_PATH),
            )

            # get_or_create_collection: if collection exists, use it;
            # if not, create it fresh — safe to call on every startup
            self.collection = self.client.get_or_create_collection(
                name=COLLECTION_NAME,
                metadata={"description": "BioAgent bioinformatics knowledge base"}
            )

            logger.info(
                f"Connected to collection '{COLLECTION_NAME}' "
                f"({self.collection.count()} documents stored)"
            )

        except Exception as e:
            logger.error(f"Failed to initialise ChromaDB: {e}")
            raise

    def add_documents(
        self,
        documents: list[str],
        embeddings: list[list[float]],
        ids: list[str],
        metadatas: list[dict] | None = None
    ) -> None:
        """
        Store text chunks and their embeddings in ChromaDB.

        Args:
            documents: List of raw text chunks to store.
            embeddings: Corresponding embedding vectors (from BioEmbedder).
            ids: Unique string ID for each document chunk.
                 Example: ["paper_001_chunk_0", "paper_001_chunk_1"]
            metadatas: Optional list of dicts with extra info per chunk.
                       Example: [{"source": "lecture_1.pdf", "page": 3}]

        Raises:
            ValueError: If list lengths don't match.
        """
        if not (len(documents) == len(embeddings) == len(ids)):
            raise ValueError(
                f"documents ({len(documents)}), embeddings ({len(embeddings)}), "
                f"and ids ({len(ids)}) must all be the same length."
            )

        logger.info(f"Adding {len(documents)} document(s) to vector store...")

        try:
            self.collection.add(
                documents=documents,    # raw text — stored for retrieval
                embeddings=embeddings,  # vectors — used for similarity search
                ids=ids,                # unique identifiers
                metadatas=metadatas or [{}] * len(documents)
            )
            logger.info(f"Successfully stored {len(documents)} chunks.")

        except Exception as e:
            logger.error(f"Failed to add documents to ChromaDB: {e}")
            raise

    def query(
        self,
        query_embedding: list[float],
        n_results: int = 3
    ) -> list[dict]:
        """
        Find the most relevant text chunks for a given query embedding.

        ChromaDB computes cosine similarity between the query vector
        and all stored vectors, returning the closest matches.

        Args:
            query_embedding: Embedding vector of the search query.
            n_results: Number of results to return (default 3).

        Returns:
            List of dicts, each containing:
            - "text": the matching text chunk
            - "id": the chunk's unique ID
            - "distance": similarity score (lower = more similar)
            - "metadata": any stored metadata for the chunk
        """
        logger.info(f"Querying vector store for top {n_results} results...")

        try:
            results = self.collection.query(
                query_embeddings=[query_embedding],  # ChromaDB expects a list
                n_results=n_results,
                include=["documents", "distances", "metadatas"]
            )

            # Unpack ChromaDB's nested response format into clean dicts
            output = []
            for i in range(len(results["documents"][0])):
                output.append({
                    "text": results["documents"][0][i],
                    "id": results["ids"][0][i],
                    "distance": results["distances"][0][i],
                    "metadata": results["metadatas"][0][i]
                })

            return output

        except Exception as e:
            logger.error(f"Failed to query ChromaDB: {e}")
            raise

    def count(self) -> int:
        """Return total number of documents stored in the collection."""
        return self.collection.count()