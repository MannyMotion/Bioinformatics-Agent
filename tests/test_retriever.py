"""
test_retriever.py

Unit tests for the RAG retriever module.
Author:  Emmanuel Ogbu (Manny)
Date:    2026-04-23
"""

import pytest
from bioagent.rag.retriever import BioRetriever, RetrievalResult


@pytest.fixture
def retriever() -> BioRetriever:
    """Initialise retriever for tests."""
    return BioRetriever()


def test_retriever_returns_results(retriever: BioRetriever) -> None:
    """Should return results for a relevant bioinformatics question."""
    results = retriever.retrieve("What is a Phred quality score?")
    assert len(results) > 0
    assert all(isinstance(r, RetrievalResult) for r in results)


def test_retriever_relevance_score(retriever: BioRetriever) -> None:
    """Relevance score should be between 0 and 100."""
    results = retriever.retrieve("DNA sequencing quality")
    for r in results:
        assert 0 <= r.relevance <= 100


def test_retriever_context_string(retriever: BioRetriever) -> None:
    """retrieve_as_context should return a non-empty string."""
    context = retriever.retrieve_as_context("What is RNA sequencing?")
    assert isinstance(context, str)
    assert len(context) > 0


def test_retriever_empty_question(retriever: BioRetriever) -> None:
    """Empty question should raise ValueError."""
    with pytest.raises(ValueError):
        retriever.retrieve("")