"""
Utility Functions for Jekyll Blog Build System

This module provides common utility functions used across the build system,
primarily for logging configuration and shared helper functions.
"""

import os
import logging
import sys
from typing import Optional, List, Tuple

logger = logging.getLogger(__name__)
    

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
  
  
def find_code_block_pairs(lines: List[str]) -> List[Tuple[int, int]]:
    """
    Find pairs of ``` lines that form code blocks.
    
    Returns:
        List of (start_line, end_line) tuples
    """
    tick_lines = []
    for i, line in enumerate(lines):
        if line.strip().startswith('```'):
            tick_lines.append(i)
    
    if len(tick_lines) % 2 != 0:
        return []  # Can't pair odd number of ticks
    
    return [(tick_lines[i], tick_lines[i + 1]) for i in range(0, len(tick_lines), 2)]


def add_output_labels(lines: List[str], label: str = "output") -> List[str]:
    """
    Add 'output' labels to code blocks that have no language specification.
    
    Args:
        lines: List of lines from the markdown file
        label: label to apply to markdown blocks which lack a language
        
    Returns:
        List of modified lines
    """
    
    # Make a copy to avoid modifying the original
    modified_lines = lines.copy()
    
    # Find all code block pairs
    code_blocks = find_code_block_pairs(lines)
    
    modifications_made = 0
    
    for start_line, end_line in code_blocks:
        opening_line = lines[start_line].strip()
        
        # Check if it's exactly "```" with no language
        if opening_line == "```":
            # Replace with "```output"
            # Preserve any leading whitespace from the original line
            original_line = lines[start_line]
            leading_whitespace = original_line[:len(original_line) - len(original_line.lstrip())]
            modified_lines[start_line] = leading_whitespace + f"```{label}\n"
            modifications_made += 1
            logger.debug(f"Added 'output' label to code block at line {start_line + 1}")
    
    if modifications_made > 0:
        logger.info(f"Added 'output' labels to {modifications_made} code block(s)")
    else:
        logger.debug("No unlabeled code blocks found")
    
    return modified_lines


def patch_markdown(markdown_path: str) -> None:
    """
    Patch a markdown file to add 'output' labels to unlabeled code blocks.
    
    Args:
        markdown_path: Path to the markdown file to patch
    """
    logger = logging.getLogger(__name__)
    
    if not os.path.exists(markdown_path):
        logger.warning(f"Markdown file not found: {markdown_path}")
        return
    
    try:
        # Read the file
        with open(markdown_path, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        # Apply the output labels
        modified_lines = add_output_labels(lines)
        
        # Write back to file if modifications were made
        if modified_lines != lines:
            with open(markdown_path, 'w', encoding='utf-8') as f:
                f.writelines(modified_lines)
            logger.debug(f"Patched markdown file: {markdown_path}")
        else:
            logger.debug(f"No changes needed for: {markdown_path}")
            
    except Exception as e:
        logger.error(f"Error patching markdown file {markdown_path}: {e}")
        raise
