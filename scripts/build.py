#!/usr/bin/env python3
"""
Jekyll Blog Build System

This script builds Jekyll blog posts from source files in multiple formats:
- R Markdown (.Rmd) files using knitr
- Quarto (.qmd) files using Quarto with knitr/jupyter engines

The script processes files in the _source directory that have proper
YYYY-MM-DD date prefixes and outputs them to the _posts directory.
"""

import os
import glob
import subprocess
import sys
import re
import logging
from pathlib import Path
from typing import List, Tuple

from utils import setup_logging


def check_working_directory() -> None:
    """
    Ensure we're running from the shackett.github.io directory.
    
    Raises:
        SystemExit: If not in the correct directory
    """
    logger = logging.getLogger(__name__)
    
    cwd = Path.cwd()
    expected_markers = ['_config.yml', '_posts', '_source']
    
    # Check if we're in the right directory
    if cwd.name != 'shackett.github.io':
        logger.warning(f"Expected to be in 'shackett.github.io' directory, but currently in '{cwd.name}'")
    
    # Check for expected Jekyll/blog structure
    missing_markers = [marker for marker in expected_markers if not (cwd / marker).exists()]
    
    if missing_markers:
        logger.error(f"Missing expected files/directories: {missing_markers}")
        logger.error("Please run this script from the shackett.github.io root directory")
        sys.exit(1)
    
    logger.debug(f"Working directory confirmed: {cwd}")


def has_date_prefix(filename: str) -> bool:
    """
    Check if filename starts with YYYY-MM-DD pattern.
    
    Args:
        filename: Filename to check
        
    Returns:
        bool: True if filename has proper date prefix
    """
    basename = os.path.basename(filename)
    # Match YYYY-MM-DD at start of filename
    date_pattern = r'^\d{4}-\d{2}-\d{2}-'
    return bool(re.match(date_pattern, basename))


def get_files_to_process() -> Tuple[List[str], List[str]]:
    """
    Get list of source files to process, filtering by date prefix.
    
    Returns:
        Tuple[List[str], List[str]]: (files_to_process, skipped_files)
    """
    logger = logging.getLogger(__name__)
    
    # Find all Rmd and qmd files in _source (top level only)
    source_files = []
    
    # Get .Rmd files
    rmd_files = glob.glob("_source/*.Rmd")
    source_files.extend(rmd_files)
    
    # Get .qmd files  
    qmd_files = glob.glob("_source/*.qmd")
    source_files.extend(qmd_files)
    
    logger.debug(f"Found {len(source_files)} total source files")
    
    # Separate files with and without date prefixes
    files_to_process = []
    skipped_files = []
    
    for file in source_files:
        if has_date_prefix(file):
            files_to_process.append(file)
        else:
            skipped_files.append(file)
    
    return files_to_process, skipped_files


def build_one(input_file: str, output_file: str, verbose: bool = False) -> None:
    """
    Build a single file using the appropriate renderer.
    
    Args:
        input_file: Path to source file (.Rmd or .qmd)
        output_file: Path where output .md file should be saved
        verbose: If True, show subprocess output; if False, capture it
        
    Raises:
        subprocess.CalledProcessError: If build process fails
    """
    logger = logging.getLogger(__name__)
    
    logger.info(f"Building {input_file}")
    
    # Get file extension to determine build method
    file_ext = Path(input_file).suffix.lower()
    
    try:
        if file_ext == ".rmd":
            # Use R build process
            subprocess.run([
                "Rscript", "scripts/build_one_Rmd.R", input_file, output_file
            ], check=True, capture_output=not verbose, text=True)
            
        elif file_ext == ".qmd":
            # Use Quarto build process
            subprocess.run([
                "python3", "scripts/build_one_quarto.py", input_file, output_file
            ], check=True, capture_output=not verbose, text=True)
            
        else:
            logger.warning(f"Unknown file type {file_ext} for {input_file}")
            return
            
        logger.debug(f"Successfully built {input_file}")
            
    except subprocess.CalledProcessError as e:
        # Clean up failed output file
        if os.path.exists(output_file):
            os.unlink(output_file)
        
        logger.error(f"Failed to compile {input_file} to {output_file}")
        logger.error(f"Error: {e}")
        if e.stdout:
            logger.debug(f"Stdout: {e.stdout}")
        if e.stderr:
            logger.error(f"Stderr: {e.stderr}")
        sys.exit(1)


def main() -> None:
    """
    Main build function.
    
    Processes all source files with proper date prefixes and builds them
    into Jekyll-compatible markdown in the _posts directory.
    """
    logger = logging.getLogger(__name__)
    
    # Check we're in the right directory
    check_working_directory()
    
    logger.info("Force rebuild enabled - regenerating all files")
    
    # Get files to process and skipped files
    files_to_process, skipped_files = get_files_to_process()
    
    # Report warnings for skipped files
    if skipped_files:
        logger.warning("The following files don't have YYYY-MM-DD date prefixes and will be skipped:")
        for file in skipped_files:
            basename = os.path.basename(file)
            logger.warning(f"  - {file}")
            logger.info(f"    Consider renaming to: YYYY-MM-DD-{basename}")
    
    if not files_to_process:
        if skipped_files:
            logger.info("No files with proper date prefixes found to build.")
        else:
            logger.info("No source files found to build.")
        return
    
    logger.info(f"Found {len(files_to_process)} source file(s) to process.")
    
    # Build each file
    for input_file in files_to_process:
        # Create output filename in _posts directory
        input_path = Path(input_file)
        output_filename = input_path.stem + ".md"
        output_file = os.path.join("_posts", output_filename)
        
        # Ensure _posts directory exists
        os.makedirs("_posts", exist_ok=True)
        
        # Build the file
        build_one(input_file, output_file)
    
    logger.info("Build completed successfully!")


if __name__ == "__main__":
    setup_logging()
    main()
