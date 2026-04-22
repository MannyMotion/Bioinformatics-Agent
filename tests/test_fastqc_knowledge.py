"""Quick test to verify FastQC lecture notes are queryable."""

from bioagent.rag.retriever import BioRetriever

retriever = BioRetriever()

questions = [
    "What does FastQC measure?",
    "What is Kmer content in sequencing?",
    "How do you interpret per base sequence quality?"
]

for question in questions:
    print(f"\nQ: {question}")
    context = retriever.retrieve_as_context(question, n_results=1)
    print(context[:300])
    print("-" * 50)