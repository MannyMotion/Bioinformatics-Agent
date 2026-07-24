"""
test_api.py

Tests the FastAPI backend endpoints for the two-phase agentic flow:
    Phase 1  POST /upload         -> detect + profile + propose questions
    Phase 2  POST /analyze/{id}   -> run the pipeline for a question

Run with the server already running on port 8000:
    python -m uvicorn bioagent.api.main:app --port 8000

Note: /analyze invokes the local Ollama model, so each test can take
30-60 seconds.

Author: Emmanuel Ogbu (Manny)
Date:   2026-07-02
"""

import requests

API = "http://localhost:8000"


def test_health():
    """Health endpoint returns ok."""
    r = requests.get(f"{API}/health")
    assert r.status_code == 200
    assert r.json()["status"] == "ok"


def test_fasta_flow():
    """Phase 1 profiles a FASTA file; Phase 2 runs the QC pipeline."""
    with open("data/sample/test.fasta", "rb") as f:
        up = requests.post(f"{API}/upload", files={"file": ("test.fasta", f)})
    assert up.status_code == 200
    data = up.json()
    assert data["file_type"] == "FASTA"
    assert data["profile"]["suggested_questions"], "agent should propose questions"

    an = requests.post(f"{API}/analyze/{data['job_id']}",
                       data={"question": ""}, timeout=600)  # "" = agent decides
    assert an.status_code == 200
    result = an.json()
    assert result["pipeline"] == "FASTA QC Pipeline"
    assert "stats" in result
    assert result["question"], "the answered question should be recorded"


def test_vcf_flow():
    """Phase 1 profiles a VCF; Phase 2 runs the variant pipeline."""
    with open("data/sample/test.vcf", "rb") as f:
        up = requests.post(f"{API}/upload", files={"file": ("test.vcf", f)})
    assert up.status_code == 200
    job_id = up.json()["job_id"]

    an = requests.post(f"{API}/analyze/{job_id}",
                       data={"question": "Are there pathogenic variants?"}, timeout=600)
    assert an.status_code == 200
    result = an.json()
    assert result["pipeline"] == "Variant Annotation Pipeline"
    assert result["stats"]["pathogenic_count"] >= 0


def test_csv_group_inference():
    """Phase 1 auto-infers control/treatment groups from column names."""
    with open("data/sample/counts.csv", "rb") as f:
        up = requests.post(f"{API}/upload", files={"file": ("counts.csv", f)})
    assert up.status_code == 200
    profile = up.json()["profile"]
    # counts.csv columns are healthy_* / cancer_* — should map to two groups.
    assert len(profile["groups"]) == 2
    assert not profile["needs_group_confirmation"]


if __name__ == "__main__":
    test_health()
    print("PASS: health")
    test_csv_group_inference()
    print("PASS: CSV group inference")
    test_fasta_flow()
    print("PASS: FASTA flow")
    test_vcf_flow()
    print("PASS: VCF flow")
    print("All API tests passed.")
