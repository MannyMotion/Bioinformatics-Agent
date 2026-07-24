"""
job_store.py

Purpose: SQLite-backed persistence for completed analysis jobs.
         Replaces the previous in-memory dict so that:
         - Jobs survive a server restart (the Q&A box keeps working).
         - Memory does not grow without bound over a long-running server.

         Each completed /upload stores its packaged JSON response here,
         keyed by job_id. The /ask endpoint reads it back to answer
         follow-up questions.

         Uses the sqlite3 standard library (no extra dependency). A fresh
         connection is opened per operation so it is safe to call from the
         FastAPI event loop and from pipeline worker threads alike.

Author:  Emmanuel Ogbu (Manny)
Date:    2026-07-02
"""

import json
import sqlite3
import time
from pathlib import Path

from bioagent.utils.logger import get_logger

logger = get_logger(__name__)

# Database file — lives at the project root alongside uploads/ and outputs/.
DB_PATH = Path("jobs.db")


def init_db() -> None:
    """Create the jobs table if it does not already exist."""
    with sqlite3.connect(DB_PATH) as conn:
        conn.execute(
            """
            CREATE TABLE IF NOT EXISTS jobs (
                job_id     TEXT PRIMARY KEY,
                data       TEXT NOT NULL,
                created_at REAL NOT NULL
            )
            """
        )
        conn.commit()
    logger.info(f"Job store ready at {DB_PATH}.")


def save_job(job_id: str, data: dict) -> None:
    """
    Persist a completed job's result payload.

    Args:
        job_id: Unique job identifier.
        data: JSON-serialisable dict (the packaged /upload response).
    """
    with sqlite3.connect(DB_PATH) as conn:
        conn.execute(
            "INSERT OR REPLACE INTO jobs (job_id, data, created_at) "
            "VALUES (?, ?, ?)",
            (job_id, json.dumps(data), time.time()),
        )
        conn.commit()


def get_job(job_id: str) -> dict | None:
    """
    Retrieve a stored job's result payload.

    Args:
        job_id: The job ID to look up.

    Returns:
        The stored dict, or None if the job_id is unknown.
    """
    with sqlite3.connect(DB_PATH) as conn:
        row = conn.execute(
            "SELECT data FROM jobs WHERE job_id = ?", (job_id,)
        ).fetchone()

    if row is None:
        return None
    return json.loads(row[0])
