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
from bioagent.agent.detector import detect_file_type
from bioagent.agent.analyst import (
    profile_file, plan_analysis, autopick_question, DataProfile
)
from bioagent.agent.explainer import explain_results
from bioagent.api import job_store
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
STATIC_DIR = Path("static")
UPLOAD_DIR.mkdir(exist_ok=True)
OUTPUT_DIR.mkdir(exist_ok=True)
STATIC_DIR.mkdir(exist_ok=True)

# Persistent job store — results are saved to SQLite (see job_store.py)
# so they survive a restart and memory does not grow without bound.
job_store.init_db()

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
# Serve vendored JS libraries (Chart.js) so plots render fully offline
app.mount("/static", StaticFiles(directory="static"), name="static")


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
    file: UploadFile = File(...)
) -> JSONResponse:
    """
    PHASE 1 of the agentic flow — upload, detect, and PROFILE a file.

    This does NOT run the heavy analysis. Instead the agent reads the data,
    describes what it is, infers sample groups (for expression matrices), and
    proposes biological questions the data can answer. The frontend then lets
    the user pick a suggested question, type their own, or let the agent
    decide — and calls POST /analyze/{job_id} to actually run it.

    Security checks applied:
    - Rate limited to 10 uploads per minute per IP
    - File size capped at 50MB
    - Only known bioinformatics extensions accepted

    Args:
        request: FastAPI request object (required by slowapi).
        file: The uploaded bioinformatics file.

    Returns:
        JSON with job_id, detected type, data profile, and suggested questions.
    """
    job_id = str(uuid.uuid4())[:8]
    logger.info(f"New upload job: {job_id} — file: {file.filename}")

    # Security check 0 — strip any directory components from the filename.
    # Prevents path traversal (e.g. "../../etc/x") from writing outside UPLOAD_DIR.
    # Path().name handles both "/" and "\\" separators regardless of host OS.
    safe_filename = Path((file.filename or "upload").replace("\\", "/")).name or "upload"

    # Security check 1 — validate file extension before saving
    file_extension = Path(safe_filename).suffix.lower()
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
    upload_path = UPLOAD_DIR / f"{job_id}_{safe_filename}"
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

    # Detect + profile in the thread pool (pandas/file reads are blocking).
    try:
        loop = asyncio.get_running_loop()
        detection, profile = await loop.run_in_executor(
            executor,
            lambda: _profile_upload(upload_path)
        )
        logger.info(f"Job {job_id} profiled: {detection.file_type}")
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.error(f"Profiling failed for job {job_id}: {e}")
        raise HTTPException(status_code=500, detail=f"Could not read file: {e}")

    record = {
        "job_id": job_id,
        "file_name": safe_filename,
        "file_path": str(upload_path),
        "file_type": detection.file_type,
        "confidence": detection.confidence,
        "detection_explanation": detection.explanation,
        "profile": profile.as_dict(),
        "status": "profiled",
    }
    job_store.save_job(job_id, record)
    return JSONResponse(content=record)


@app.post("/analyze/{job_id}")
@limiter.limit("10/minute")
async def analyze_job(
    request: Request,
    job_id: str,
    question: str = Form(default=""),
    control_label: str = Form(default=""),
    treatment_label: str = Form(default=""),
    control_cols: str = Form(default=""),
    treatment_cols: str = Form(default="")
) -> JSONResponse:
    """
    PHASE 2 of the agentic flow — answer a question by running the pipeline.

    Uses the profile stored at upload. If `question` is blank the agent picks
    its own best question. For expression data the caller may override the
    inferred groups by passing comma-separated `control_cols`/`treatment_cols`.

    Args:
        request: FastAPI request object (required by slowapi).
        job_id: The job ID from POST /upload.
        question: Chosen or typed question (blank = let the agent decide).
        control_label / treatment_label: Optional group display names.
        control_cols / treatment_cols: Optional explicit group columns (CSV).

    Returns:
        JSON with the full analysis results and interpretation.
    """
    record = job_store.get_job(job_id)
    if record is None:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found.")

    file_path = record.get("file_path")
    if not file_path or not Path(file_path).exists():
        raise HTTPException(
            status_code=410,
            detail="Uploaded file is no longer available. Please re-upload."
        )

    profile = DataProfile.from_dict(record.get("profile", {}))
    chosen_question = question.strip() or autopick_question(profile)

    # Optional group override from the UI (comma-separated column names).
    ctrl = [c.strip() for c in control_cols.split(",") if c.strip()]
    treat = [c.strip() for c in treatment_cols.split(",") if c.strip()]
    groups = None
    if ctrl and treat:
        cl = control_label.strip() or "control"
        tl = treatment_label.strip() or "treatment"
        groups = {cl: ctrl, tl: treat}

    plan = plan_analysis(chosen_question, profile, groups)
    logger.info(f"Job {job_id} analysing — question: {chosen_question[:60]}")

    try:
        loop = asyncio.get_running_loop()
        result, decision = await loop.run_in_executor(
            executor,
            lambda: _run_analysis(Path(file_path), plan, chosen_question)
        )
        logger.info(f"Job {job_id} complete: {decision.pipeline_name}")
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.error(f"Analysis failed for job {job_id}: {e}")
        raise HTTPException(status_code=500, detail=f"Analysis failed: {e}")

    response_data = _package_results(job_id, result, decision)
    response_data["question"] = chosen_question
    response_data["plan_steps"] = plan.steps
    response_data["profile"] = record.get("profile", {})
    response_data["file_path"] = file_path
    response_data["status"] = "complete"
    job_store.save_job(job_id, response_data)
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
    job_data = job_store.get_job(job_id)
    if job_data is None:
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


def _profile_upload(upload_path: Path) -> tuple[Any, Any]:
    """
    Detect and profile an upload (executed in a worker thread).

    Args:
        upload_path: Path to the saved upload.

    Returns:
        Tuple of (DetectionResult, DataProfile).

    Raises:
        ValueError: If the file type cannot be determined.
    """
    detection = detect_file_type(upload_path)
    if detection.file_type.upper() == "UNKNOWN":
        raise ValueError(
            f"Could not determine file type for {upload_path.name}. "
            f"Supported: FASTA, FASTQ, VCF, CSV, TSV."
        )
    profile = profile_file(upload_path, detection)
    return detection, profile


def _run_analysis(upload_path: Path, plan: Any, question: str) -> tuple[Any, Any]:
    """
    Run the full analysis for a planned question (executed in a worker thread).

    Routes the file to the pipeline named in the plan (passing any explicit
    sample-group columns), then generates a question-aware Ollama
    interpretation. The pipeline's own template interpretation is the fallback
    if Ollama is unavailable.

    Args:
        upload_path: Path to the saved upload.
        plan: The AnalysisPlan from the analyst.
        question: The question being answered (shapes the interpretation).

    Returns:
        Tuple of (pipeline_result, RoutingDecision).
    """
    params = plan.params
    result, decision = route_file(
        upload_path,
        output_dir=OUTPUT_DIR,
        use_rag=False,  # pipeline fills result.interpretation with the template
        rnaseq_control=params.get("control_label", "control"),
        rnaseq_treatment=params.get("treatment_label", "treatment"),
        rnaseq_control_cols=params.get("control_cols") or None,
        rnaseq_treatment_cols=params.get("treatment_cols") or None,
    )

    # Upgrade the template interpretation to a question-aware Ollama one,
    # keeping the template as the fallback if Ollama is unavailable.
    result.interpretation = explain_results(
        pipeline_result=result,
        pipeline_name=decision.pipeline_name,
        stats=_extract_stats(result),
        warnings=result.warnings,
        use_rag=True,
        fallback=result.interpretation,
        question=question,
    )

    return result, decision


def _extract_stats(result: Any) -> dict:
    """
    Pull the pipeline-specific summary statistics off a result object.

    Each pipeline returns a different result type; we detect which by the
    attributes present and return a flat dict suitable for JSON and for
    the Ollama prompt.

    Args:
        result: A QCResult, RNAseqResult, or VariantResult.

    Returns:
        Dictionary of key statistics (empty if the type is unrecognised).
    """
    if hasattr(result, "mean_gc"):
        return {
            "total_sequences": result.total_sequences,
            "total_bases": result.total_bases,
            "mean_gc": result.mean_gc,
            "median_length": result.median_length,
            "min_length": result.min_length,
            "max_length": result.max_length,
            "low_complexity_count": result.low_complexity_count,
        }
    elif hasattr(result, "upregulated_count"):
        return {
            "total_genes": result.total_genes,
            "genes_tested": result.genes_tested,
            "upregulated_count": result.upregulated_count,
            "downregulated_count": result.downregulated_count,
        }
    elif hasattr(result, "pathogenic_count"):
        return {
            "total_variants": result.total_variants,
            "pass_filter_count": result.pass_filter_count,
            "annotated_count": result.annotated_count,
            "pathogenic_count": result.pathogenic_count,
        }
    return {}


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
    stats = _extract_stats(result)
    if stats:
        base["stats"] = stats

    return base


    #& "C:/Users/Invate/anaconda3/envs/bioagent/python.exe" -m uvicorn bioagent.api.main:app --host 127.0.0.1 --port 8000
