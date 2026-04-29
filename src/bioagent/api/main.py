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

Inputs:  Multipart file upload via HTTP POST
Outputs: JSON responses with analysis results and plot paths

Author:  Emmanuel Ogbu (Manny)
Date:    2026-04-28
"""

import uuid
import shutil
from pathlib import Path
from typing import Any

from fastapi import FastAPI, UploadFile, File, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from fastapi.responses import JSONResponse

from bioagent.agent.router import route_file
from bioagent.utils.logger import get_logger

logger = get_logger(__name__)

# Directories for uploaded files and outputs
UPLOAD_DIR = Path("uploads")
OUTPUT_DIR = Path("outputs")
UPLOAD_DIR.mkdir(exist_ok=True)
OUTPUT_DIR.mkdir(exist_ok=True)

# In-memory job store — stores results by job_id
# In production this would be a database (SQLite/PostgreSQL)
job_store: dict[str, dict] = {}

# Initialise FastAPI app
app = FastAPI(
    title="BioAgent API",
    description="Agentic Bioinformatics Analysis System",
    version="0.3.0"
)

# Allow frontend to talk to backend (CORS)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # in production, restrict to your domain
    allow_methods=["*"],
    allow_headers=["*"],
)

# Serve output plots as static files
app.mount("/outputs", StaticFiles(directory="outputs"), name="outputs")


@app.get("/health")
def health_check() -> dict:
    """
    Health check endpoint.
    Returns server status — used by frontend to verify API is running.
    """
    return {"status": "ok", "version": "0.3.0"}


@app.post("/upload")
from fastapi import Form

async def upload_file(
    file: UploadFile = File(...),
    control_label: str = Form(default="control"),
    treatment_label: str = Form(default="treatment")
) -> JSONResponse:
    """
    Upload a bioinformatics file and run the appropriate pipeline.

    The system will:
    1. Save the uploaded file
    2. Auto-detect its type
    3. Route to the correct pipeline
    4. Run the full analysis
    5. Return results with a job_id for retrieval

    Args:
        file: The uploaded bioinformatics file (FASTA, VCF, CSV etc.)
        control_label: For CSV files — column prefix for control samples
        treatment_label: For CSV files — column prefix for treatment samples

    Returns:
        JSON with job_id, detected file type, pipeline used, and results.
    """
    # Generate unique job ID for this analysis
    job_id = str(uuid.uuid4())[:8]
    logger.info(f"New upload job: {job_id} — file: {file.filename}")

    # Save uploaded file to disk
    upload_path = UPLOAD_DIR / f"{job_id}_{file.filename}"
    try:
        with open(upload_path, "wb") as buffer:
            shutil.copyfileobj(file.file, buffer)
        logger.info(f"Saved upload: {upload_path}")
    except Exception as e:
        logger.error(f"Failed to save upload: {e}")
        raise HTTPException(status_code=500, detail=f"Failed to save file: {e}")

    # Run the pipeline via router
    try:
        result, decision = route_file(
            upload_path,
            output_dir=OUTPUT_DIR,
            use_rag=False,  # set True when Ollama is integrated
            rnaseq_control=control_label,
            rnaseq_treatment=treatment_label
        )
        logger.info(f"Job {job_id} complete: {decision.pipeline_name}")
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.error(f"Pipeline failed for job {job_id}: {e}")
        raise HTTPException(status_code=500, detail=f"Analysis failed: {e}")

    # Package results into JSON-serialisable format
    response_data = _package_results(job_id, result, decision)

    # Store in job store for later retrieval
    job_store[job_id] = response_data

    return JSONResponse(content=response_data)


@app.get("/analyse/{job_id}")
def get_results(job_id: str) -> JSONResponse:
    """
    Retrieve results for a previously submitted job.

    Args:
        job_id: The job ID returned by /upload

    Returns:
        JSON with full analysis results.
    """
    if job_id not in job_store:
        raise HTTPException(
            status_code=404,
            detail=f"Job {job_id} not found."
        )
    return JSONResponse(content=job_store[job_id])


@app.get("/jobs")
def list_jobs() -> JSONResponse:
    """List all completed analysis jobs."""
    jobs = [
        {
            "job_id": jid,
            "file_name": data.get("file_name"),
            "pipeline": data.get("pipeline"),
            "file_type": data.get("file_type")
        }
        for jid, data in job_store.items()
    ]
    return JSONResponse(content={"jobs": jobs})


def _package_results(job_id: str, result: Any, decision: Any) -> dict:
    """
    Convert pipeline result objects into JSON-serialisable dictionaries.

    Different pipelines return different result types (QCResult,
    RNAseqResult, VariantResult) — this function handles all of them
    and produces a consistent response format.

    Args:
        job_id: Unique job identifier.
        result: Pipeline result object.
        decision: RoutingDecision from the router.

    Returns:
        Dictionary safe to serialise as JSON.
    """
    # Convert plot paths to URL paths for frontend access
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

    # Add pipeline-specific fields
    if hasattr(result, "mean_gc"):
        # FASTA QC result
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
        # RNA-seq result
        base["stats"] = {
            "total_genes": result.total_genes,
            "genes_tested": result.genes_tested,
            "upregulated_count": result.upregulated_count,
            "downregulated_count": result.downregulated_count,
        }

    elif hasattr(result, "pathogenic_count"):
        # Variant result
        base["stats"] = {
            "total_variants": result.total_variants,
            "pass_filter_count": result.pass_filter_count,
            "annotated_count": result.annotated_count,
            "pathogenic_count": result.pathogenic_count,
        }

    return base