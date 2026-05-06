"""
main.py

Purpose: FastAPI backend for the BioAgent system.
         Provides REST API endpoints for file upload, analysis,
         and results retrieval. This is the server that the
         frontend communicates with.

         Endpoints:
         POST /upload     — upload a bioinformatics file
         GET  /analyse/{job_id} — get analysis results
         GET  /health     — check server is running
         POST /ask/{job_id} — ask Ollama a follow-up question

Inputs:  Multipart file upload via HTTP POST
Outputs: JSON responses with analysis results and plot paths

Author:  Emmanuel Ogbu (Manny)
Date:    2026-04-28
"""

import asyncio
import uuid
import shutil
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

from fastapi import FastAPI, UploadFile, File, HTTPException, Form, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from fastapi.responses import JSONResponse
from slowapi import Limiter, _rate_limit_exceeded_handler
from slowapi.util import get_remote_address
from slowapi.errors import RateLimitExceeded

from bioagent.agent.router import route_file
from bioagent.utils.logger import get_logger

logger = get_logger(__name__)

# Security constants
MAX_FILE_SIZE = 50 * 1024 * 1024  # 50MB — reject anything larger
ALLOWED_EXTENSIONS = {
    '.fasta', '.fa', '.fna', '.fastq',
    '.vcf', '.csv', '.tsv', '.txt'
}

# Directories for uploaded files and outputs
UPLOAD_DIR = Path("uploads")
OUTPUT_DIR = Path("outputs")
UPLOAD_DIR.mkdir(exist_ok=True)
OUTPUT_DIR.mkdir(exist_ok=True)

# In-memory job store — stores results by job_id
# In production this would be a database (SQLite/PostgreSQL)
job_store: dict[str, dict] = {}

# Rate limiter — prevents abuse of the API
# key_func=get_remote_address limits per IP address
limiter = Limiter(key_func=get_remote_address)

# Thread pool for running CPU-bound pipeline in background
# Prevents matplotlib/numpy from crashing the async event loop on Windows
executor = ThreadPoolExecutor(max_workers=2)

# Initialise FastAPI app
app = FastAPI(
    title="BioAgent API",
    description="Agentic Bioinformatics Analysis System",
    version="0.5.0"
)

# Wire rate limiter into FastAPI exception handling
app.state.limiter = limiter
app.add_exception_handler(RateLimitExceeded, _rate_limit_exceeded_handler)

# Allow frontend to talk to backend (CORS)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # in production, restrict to your domain
    allow_methods=["*"],
    allow_headers=["*"],
)

# Serve output plots and frontend as static files
app.mount("/outputs", StaticFiles(directory="outputs"), name="outputs")
app.mount("/frontend", StaticFiles(directory="frontend"), name="frontend")


@app.get("/health")
def health_check() -> dict:
    """
    Health check endpoint.
    Returns server status — used by frontend to verify API is running.
    """
    return {"status": "ok", "version": "0.5.0"}


@app.post("/upload")
@limiter.limit("10/minute")  # max 10 uploads per minute per IP
async def upload_file(
    request: Request,  # required by slowapi for rate limiting
    file: UploadFile = File(...),
    control_label: str = Form(default="control"),
    treatment_label: str = Form(default="treatment")
) -> JSONResponse:
    """
    Upload a bioinformatics file and run the appropriate pipeline.

    Security checks applied:
    - Rate limited to 10 uploads per minute per IP
    - File size capped at 50MB
    - Only known bioinformatics extensions accepted

    Args:
        request: FastAPI request object (required by slowapi).
        file: The uploaded bioinformatics file.
        control_label: For CSV files — column prefix for control samples.
        treatment_label: For CSV files — column prefix for treatment samples.

    Returns:
        JSON with job_id, detected file type, pipeline used, and results.
    """
    job_id = str(uuid.uuid4())[:8]
    logger.info(f"New upload job: {job_id} — file: {file.filename}")

    # Security check 1 — validate file extension before saving
    file_extension = Path(file.filename).suffix.lower()
    if file_extension not in ALLOWED_EXTENSIONS:
        logger.warning(f"Rejected upload: disallowed extension '{file_extension}'")
        raise HTTPException(
            status_code=400,
            detail=(
                f"File type '{file_extension}' not allowed. "
                f"Supported formats: {', '.join(sorted(ALLOWED_EXTENSIONS))}"
            )
        )

    # Save uploaded file to disk
    upload_path = UPLOAD_DIR / f"{job_id}_{file.filename}"
    try:
        with open(upload_path, "wb") as buffer:
            shutil.copyfileobj(file.file, buffer)
        logger.info(f"Saved upload: {upload_path}")
    except Exception as e:
        logger.error(f"Failed to save upload: {e}")
        raise HTTPException(status_code=500, detail=f"Failed to save file: {e}")

    # Security check 2 — validate file size after saving
    file_size = upload_path.stat().st_size
    if file_size > MAX_FILE_SIZE:
        upload_path.unlink()  # delete the oversized file immediately
        logger.warning(
            f"Rejected upload: file too large "
            f"({file_size / 1024 / 1024:.1f}MB > 50MB limit)"
        )
        raise HTTPException(
            status_code=413,
            detail=(
                f"File too large ({file_size / 1024 / 1024:.1f}MB). "
                f"Maximum allowed size is 50MB."
            )
        )

    # Run pipeline in thread pool — prevents crashes on Windows
    try:
        loop = asyncio.get_event_loop()
        result, decision = await loop.run_in_executor(
            executor,
            lambda: route_file(
                upload_path,
                output_dir=OUTPUT_DIR,
                use_rag=False,
                rnaseq_control=control_label,
                rnaseq_treatment=treatment_label
            )
        )
        logger.info(f"Job {job_id} complete: {decision.pipeline_name}")
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.error(f"Pipeline failed for job {job_id}: {e}")
        raise HTTPException(status_code=500, detail=f"Analysis failed: {e}")

    response_data = _package_results(job_id, result, decision)
    job_store[job_id] = response_data
    return JSONResponse(content=response_data)


@app.post("/ask/{job_id}")
@limiter.limit("20/minute")  # more generous limit for Q&A
async def ask_question(
    request: Request,
    job_id: str,
    question: str = Form(...)
) -> JSONResponse:
    """
    Answer a follow-up question about a completed analysis using Ollama.

    Args:
        request: FastAPI request object (required by slowapi).
        job_id: The job ID from a previous /upload call.
        question: The user's follow-up question.

    Returns:
        JSON with Ollama's answer.
    """
    if job_id not in job_store:
        raise HTTPException(
            status_code=404,
            detail=f"Job {job_id} not found."
        )

    # Security: cap question length to prevent prompt injection
    if len(question) > 500:
        raise HTTPException(
            status_code=400,
            detail="Question too long. Maximum 500 characters."
        )

    job_data = job_store[job_id]

    try:
        from bioagent.agent.explainer import answer_question

        answer = answer_question(
            question=question,
            pipeline_name=job_data.get("pipeline", ""),
            stats=job_data.get("stats", {}),
            interpretation=job_data.get("interpretation", "")
        )

        return JSONResponse(content={
            "job_id": job_id,
            "question": question,
            "answer": answer
        })

    except Exception as e:
        logger.error(f"Q&A failed for job {job_id}: {e}")
        raise HTTPException(status_code=500, detail=str(e))


def _package_results(job_id: str, result: Any, decision: Any) -> dict:
    """
    Convert pipeline result objects into JSON-serialisable dictionaries.

    Args:
        job_id: Unique job identifier.
        result: Pipeline result object.
        decision: RoutingDecision from the router.

    Returns:
        Dictionary safe to serialise as JSON.
    """
    # Convert Windows backslashes to forward slashes for URL paths
    plot_urls = [
        "/" + p.replace("\\", "/")
        for p in result.plot_paths
    ]

    base = {
        "job_id": job_id,
        "file_name": result.file_name,
        "file_type": decision.file_type,
        "pipeline": decision.pipeline_name,
        "confidence": decision.confidence,
        "reasoning": decision.reasoning,
        "warnings": result.warnings,
        "interpretation": result.interpretation,
        "plot_urls": plot_urls,
    }

    # Add pipeline-specific stats
    if hasattr(result, "mean_gc"):
        base["stats"] = {
            "total_sequences": result.total_sequences,
            "total_bases": result.total_bases,
            "mean_gc": result.mean_gc,
            "median_length": result.median_length,
            "min_length": result.min_length,
            "max_length": result.max_length,
            "low_complexity_count": result.low_complexity_count,
        }
    elif hasattr(result, "upregulated_count"):
        base["stats"] = {
            "total_genes": result.total_genes,
            "genes_tested": result.genes_tested,
            "upregulated_count": result.upregulated_count,
            "downregulated_count": result.downregulated_count,
        }
    elif hasattr(result, "pathogenic_count"):
        base["stats"] = {
            "total_variants": result.total_variants,
            "pass_filter_count": result.pass_filter_count,
            "annotated_count": result.annotated_count,
            "pathogenic_count": result.pathogenic_count,
        }

    return base