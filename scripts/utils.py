"""
Utility Functions for Jekyll Blog Build System

This module provides common utility functions used across the build system,
primarily for logging configuration and shared helper functions.
"""

import logging
import sys
from typing import Optional


def setup_logging(level: Optional[int] = None) -> None:
    """
    Set up logging configuration with debug level.
    
    Configures the root logger with a consistent format and debug-level
    output to stdout. This provides detailed information about the build
    process for troubleshooting.
    
    Args:
        level: Optional logging level. Defaults to DEBUG if not specified.
    """
    if level is None:
        level = logging.DEBUG
        
    logging.basicConfig(
        level=level,
        format='%(levelname)s: %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout)
        ]
    )


def get_logger(name: str) -> logging.Logger:
    """
    Get a logger instance for the specified module.
    
    Args:
        name: Usually __name__ from the calling module
        
    Returns:
        logging.Logger: Configured logger instance
    """
    return logging.getLogger(name)
