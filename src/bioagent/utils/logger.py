"""
logger.py

Purpose: Centralised logging configuration for the BioAgent system.
         All modules import get_logger() from here — ensures consistent
         log format across the entire application.

Inputs:  Module name (passed as __name__ from calling module)
Outputs: Configured logging.Logger instance

Author:  Emmanuel Ogbu (Manny)
Date:    2026-04-22
"""

import logging
import sys
from pathlib import Path
dine
LOG_FILE = Path(__file__).parent.parent.parent.parent / "bioagent.log"
LOG_FORMAT = "%(asctime)s | %(levelname)-8s | %(name)s | %(message)s"
DATE_FORMAT = "%Y-%m-%d %H:%M:%S"


def get_logger(name: str) -> logging.Logger:
    """
    Create and configure a logger for a given module.

    Args:
        name: The module name — always pass __name__ when calling.

    Returns:
        Configured Logger instance ready to use.
    """
    logger = logging.getLogger(name)

    # Guard: don't add handlers if already configured
    if logger.handlers:
        return logger

    logger.setLevel(logging.DEBUG)

    # Console handler — INFO and above
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(logging.Formatter(LOG_FORMAT, DATE_FORMAT))

    # File handler — DEBUG and above, captures everything
    file_handler = logging.FileHandler(LOG_FILE, encoding="utf-8")
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(logging.Formatter(LOG_FORMAT, DATE_FORMAT))

    logger.addHandler(console_handler)
    logger.addHandler(file_handler)

    return logger