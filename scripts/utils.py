"""
Utility Functions for Jekyll Blog Build System

This module provides common utility functions used across the build system,
primarily for logging configuration and shared helper functions.
"""

import re
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
        modified_lines = ensure_code_block_spacing(modified_lines)
        
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


def ensure_code_block_spacing(lines: List[str]) -> List[str]:
    """
    Ensure proper spacing around code blocks and remove consecutive empty lines.
    
    Args:
        lines: List of lines from the markdown file
        
    Returns:
        List of lines with proper code block spacing
    """
    logger = logging.getLogger(__name__)
    
    if not lines:
        return lines
    
    # First pass: ensure newlines around code blocks
    result = []
    
    for i, line in enumerate(lines):
        stripped = line.strip()
        
        # Check if this is a code block opener (``` with optional language)
        if stripped.startswith('```') and stripped != '```':
            # Ensure newline before opening code block (if not at start of file)
            if i > 0 and result and result[-1].strip():
                result.append('\n')
            result.append(line)
        
        # Check if this is a code block closer (just ```)
        elif stripped == '```':
            result.append(line)
            # Ensure newline after closing code block (if not at end of file)
            if i < len(lines) - 1 and lines[i + 1].strip():
                result.append('\n')
        
        else:
            result.append(line)
    
    # Second pass: remove consecutive empty lines
    final_result = []
    prev_was_empty = False
    
    for line in result:
        is_empty = not line.strip()
        
        if is_empty:
            if not prev_was_empty:
                final_result.append(line)
            prev_was_empty = True
        else:
            final_result.append(line)
            prev_was_empty = False
    
    logger.debug(f"Code block spacing: {len(lines)} -> {len(final_result)} lines")
    return final_result


def unindent_html_output(lines: List[str], start_line: int, end_line: int) -> None:
    """
    Remove indentation from HTML content within display blocks.
    
    Detects HTML content and removes the minimum indentation level from all lines
    in the block. This handles both direct display() calls (4 spaces) and nested
    calls from within functions (8+ spaces).
    
    Args:
        lines: List of all lines in the content
        start_line: Starting line index of the block content (after opener)
        end_line: Ending line index of the block content (before closer)
        
    Modifies:
        lines: Removes indentation from HTML content in-place
    """
    content_lines = lines[start_line + 1:end_line]
    
    # Check if content contains HTML tags
    contains_html = any(
        re.search(r'<(div|table|span|p|h[1-6]|style|figcaption|img|a)\b', line, re.IGNORECASE) 
        for line in content_lines if line.strip()
    )
    
    if not contains_html:
        return  # Not HTML content, leave indentation intact
    
    # Find non-empty lines to calculate minimum indentation
    non_empty_lines = [line for line in content_lines if line.strip()]
    if not non_empty_lines:
        return  # No content to process
    
    # Calculate minimum indentation level
    min_indent = min(len(line) - len(line.lstrip()) for line in non_empty_lines)
    
    if min_indent == 0:
        return  # Already unindented
    
    # Remove the detected indentation from all lines
    for i in range(start_line + 1, end_line):
        if lines[i].startswith(' ' * min_indent):
            lines[i] = lines[i][min_indent:]
        elif lines[i].strip():  # Non-empty line with less indentation
            lines[i] = lines[i].lstrip()  # Remove all leading whitespace
