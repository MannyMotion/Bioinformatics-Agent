"""
test_api.py

Tests the FastAPI backend endpoints.
Run with the server already running on port 8000.
Author: Emmanuel Ogbu (Manny)
Date:   2026-04-28
"""

import urllib.request
import json


def test_health():
    """Test health endpoint."""
    response = urllib.request.urlopen("http://localhost:8000/health")
    result = json.loads(response.read())
    print(f"Health: {result}")
    assert result["status"] == "ok"
    print("PASS: health check")


def test_upload_fasta():
    """Test uploading a FASTA file."""
    with open("data/sample/test.fasta", "rb") as f:
        file_data = f.read()

    boundary = "boundary123"
    body = (
        b"--" + boundary.encode() + b"\r\n"
        b"Content-Disposition: form-data; name=\"file\"; filename=\"test.fasta\"\r\n"
        b"Content-Type: text/plain\r\n\r\n"
        + file_data + b"\r\n"
        b"--" + boundary.encode() + b"--\r\n"
    )

    req = urllib.request.Request(
        "http://localhost:8000/upload",
        data=body,
        headers={"Content-Type": f"multipart/form-data; boundary={boundary}"}
    )

    response = urllib.request.urlopen(req)
    result = json.loads(response.read())

    print(f"Job ID: {result['job_id']}")
    print(f"Pipeline: {result['pipeline']}")
    print(f"File type: {result['file_type']}")
    print(f"Stats: {result['stats']}")
    print(f"Plots: {result['plot_urls']}")
    assert result["pipeline"] == "FASTA QC Pipeline"
    print("PASS: FASTA upload")


def test_upload_vcf():
    """Test uploading a VCF file."""
    with open("data/sample/test.vcf", "rb") as f:
        file_data = f.read()

    boundary = "boundary456"
    body = (
        b"--" + boundary.encode() + b"\r\n"
        b"Content-Disposition: form-data; name=\"file\"; filename=\"test.vcf\"\r\n"
        b"Content-Type: text/plain\r\n\r\n"
        + file_data + b"\r\n"
        b"--" + boundary.encode() + b"--\r\n"
    )

    req = urllib.request.Request(
        "http://localhost:8000/upload",
        data=body,
        headers={"Content-Type": f"multipart/form-data; boundary={boundary}"}
    )

    response = urllib.request.urlopen(req)
    result = json.loads(response.read())

    print(f"Job ID: {result['job_id']}")
    print(f"Pipeline: {result['pipeline']}")
    print(f"Pathogenic variants: {result['stats']['pathogenic_count']}")
    assert result["pipeline"] == "Variant Annotation Pipeline"
    print("PASS: VCF upload")


if __name__ == "__main__":
    test_health()
    print()
    test_upload_fasta()
    print()
    test_upload_vcf()
    print()
    print("All API tests passed.")