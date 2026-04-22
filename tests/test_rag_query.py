"""Quick test to verify RAG retrieval is working."""

from bioagent.rag.embedder import BioEmbedder
from bioagent.rag.vector_store import BioVectorStore

embedder = BioEmbedder()
store = BioVectorStore()

query = "What is a Phred quality score?"
query_vector = embedder.embed_single(query)
results = store.query(query_vector, n_results=2)

print(f"Query: {query}")
print()
for i, r in enumerate(results):
    distance = r["distance"]
    text = r["text"][:200]
    print(f"Result {i+1} (distance: {distance:.4f}):")
    print(text)
    print()