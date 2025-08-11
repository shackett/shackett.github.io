"""
Quarto Markdown Cleaning Module

This module contains functions to clean and process Quarto-generated markdown
to make it compatible with Jekyll.
"""

from typing import List, Tuple

import re
import logging
from pathlib import Path

logger = logging.getLogger(__name__)

def clean_jekyll_output(file_path: str) -> None:
    """
    Main function to clean up Quarto output for Jekyll compatibility.
    
    This function orchestrates all the cleaning steps:
    1. Remove Quarto artifacts
    2. Process tables
    3. Process figures
    4. Process code blocks
    5. Clean up spacing
    
    Args:
        file_path: Path to the markdown file to clean
        
    Raises:
        FileNotFoundError: If the input file doesn't exist
        PermissionError: If unable to write to the file
    """
    
    with open(file_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    filename_stem = Path(file_path).stem
    
    # Process in stages
    logger.debug("Processing Quarto artifacts...")
    content = remove_quarto_artifacts(content)
    
    logger.debug("Processing tables...")
    content = process_quarto_table(content)
    
    logger.debug("Processing figures...")
    content = process_quarto_figures(content, filename_stem)
    
    logger.debug("Processing code blocks...")
    content = process_code_blocks_and_output(content)
    
    # Remove trailing empty lines
    lines = content.split('\n')
    while lines and not lines[-1].strip():
        lines.pop()
    content = '\n'.join(lines)
    
    with open(file_path, 'w', encoding='utf-8') as f:
        f.write(content)
    
    logger.debug(f"Cleaned Jekyll output: {file_path}")
    logger.debug(f"Cleaned content preview:\n{content[:500]}")

def remove_quarto_artifacts(content: str) -> str:
    """
    Remove basic Quarto cell artifacts including Jupyter-specific headers.
    
    Args:
        content: Raw markdown content from Quarto
        
    Returns:
        str: Content with Quarto artifacts removed
    """
    
    # Convert code block headers to simple format first
    # ``` {.python .cell-code} -> ```python
    content = re.sub(r'``` \{\.python[^}]*\}', '```python', content)
    content = re.sub(r'``` \{\.r[^}]*\}', '```r', content)
    
    # Remove ALL Quarto cell divider blocks (general approach)
    # This catches any ::: {...} or :::: {...} blocks that contain cell-related classes
    # Pattern: one or more colons, optional whitespace, {anything with "cell" in it}, optional whitespace, newline
    content = re.sub(r'^:::+\s*\{[^}]*\.cell[^}]*\}\s*$', '', content, flags=re.MULTILINE)
    
    # Remove general cell output blocks
    # Pattern: ::: {anything with "cell-output"} 
    content = re.sub(r'^:::+\s*\{[^}]*cell-output[^}]*\}\s*$', '', content, flags=re.MULTILINE)
    
    # Remove code block output wrappers
    # ``` {.cell-output...} 
    content = re.sub(r'^``` \{[^}]*cell-output[^}]*\}\s*$', '', content, flags=re.MULTILINE)
    
    # Remove simple cell dividers (::::: cell, etc.)
    content = re.sub(r'^:::+\s*cell\s*$', '', content, flags=re.MULTILINE)
    
    # Clean up orphaned ::: markers (any standalone colons that might be left)
    content = re.sub(r'^:::+\s*$', '', content, flags=re.MULTILINE)
    
    # Remove excessive blank lines that might result from the cleanup
    content = re.sub(r'\n\s*\n\s*\n', '\n\n', content)
    
    return content


def process_quarto_table(content: str) -> str:
    """
    Process Quarto table output for Jekyll compatibility.
    
    Handles different types of table output from Quarto:
    - R kable tables wrapped in cell-output-display blocks
    - Jupyter HTML display tables  
    - Table caption formatting issues
    
    Args:
        content: Markdown content containing Quarto table output
        
    Returns:
        str: Content with properly formatted tables
    """
    
    # Process different table types in order
    content = _process_kable_tables(content)            # Handle R kable in cell blocks
    content = _process_table_captions(content)          # Fix caption formatting
    content = _process_jupyter_display_tables(content)  # Cleanup: wrap unwrapped tables
    content = _remove_table_stylers(content)            # remove styler CSS
    content = _align_text_tables(content)               # Fix alignment of markdown titles
    
    return content
  

def process_quarto_figures(content: str, filename_stem: str) -> str:
    """
    Process figure paths and display blocks.
    
    Converts Quarto figure paths to Jekyll-compatible paths and removes
    Quarto-specific display wrappers.
    
    Args:
        content: Markdown content containing figure references
        filename_stem: Base filename for constructing figure paths
        
    Returns:
        str: Content with corrected figure paths
    """
    
    # Fix figure paths - convert Quarto's figure-markdown paths to our expected paths
    content = re.sub(
        rf'{re.escape(filename_stem)}_files/figure-markdown/',
        f'figure/source/{filename_stem}/',
        content
    )
    
    # Handle figure display blocks - remove the wrapper but keep the image
    # Pattern: ::: cell-output-display\n![caption](path)\n:::
    def replace_figure_display(match):
        figure_content = match.group(1).strip()
        # If it's an image, keep just the image
        if figure_content.startswith('!['):
            return figure_content
        return figure_content
    
    pattern = r'::: cell-output-display\s*\n(!\[.*?\].*?)\s*\n:::'
    content = re.sub(pattern, replace_figure_display, content, flags=re.DOTALL)
    
    # Also handle cases where cell-output-display appears as a standalone line
    content = re.sub(r'\n::: cell-output-display\s*\n', '\n', content)
    content = re.sub(r'``` cell-output-display\s*\n', '', content)
    content = re.sub(r'``` \{\.cell-output-display\}\s*\n', '', content)
    
    # Convert relative figure paths to Jekyll-absolute paths
    # Jekyll interprets figure/source/... as relative to post URL, but we want site root
    content = re.sub(r'!\[([^\]]*)\]\(figure/source/', r'![\1](/figure/source/', content)
    
    return content


def process_code_blocks_and_output(content: str) -> str:
    """
    Process code blocks and their output with proper spacing.
    
    Ensures proper spacing between code blocks and their output,
    and formats output in Jekyll-compatible code blocks.
    
    Args:
        content: Markdown content containing code blocks
        
    Returns:
        str: Content with properly formatted code blocks
    """
    lines = content.split('\n')
    cleaned_lines = []
    i = 0
    
    while i < len(lines):
        line = lines[i]
        
        # Handle code blocks
        if line.strip() == '```python' or line.strip() == '```r':
            # Add the opening code block
            cleaned_lines.append(line)
            i += 1
            
            # Add all lines until closing ```
            while i < len(lines) and lines[i].strip() != '```':
                cleaned_lines.append(lines[i])
                i += 1
            
            # Add closing ```
            if i < len(lines):
                cleaned_lines.append(lines[i])
                i += 1
            
            # Check if next non-empty line is output (starts with spaces)
            next_line_idx = i
            while next_line_idx < len(lines) and not lines[next_line_idx].strip():
                next_line_idx += 1
            
            if (next_line_idx < len(lines) and 
                lines[next_line_idx].startswith('    ') and 
                not lines[next_line_idx].strip().startswith('```')):
                
                # Add blank line before output
                cleaned_lines.append('')
                
                # Convert output to proper Jekyll format
                cleaned_lines.append('```')  # Start output block
                
                # Add the output line(s) without leading spaces
                while (next_line_idx < len(lines) and 
                       lines[next_line_idx].startswith('    ')):
                    # Remove leading 4 spaces from output
                    output_line = lines[next_line_idx][4:] if len(lines[next_line_idx]) > 4 else lines[next_line_idx].strip()
                    cleaned_lines.append(output_line)
                    next_line_idx += 1
                
                cleaned_lines.append('```')  # End output block
                i = next_line_idx
                
                # Add blank line after output if there's more content
                if i < len(lines):
                    cleaned_lines.append('')
            else:
                # Add blank line after code block if no output follows
                if i < len(lines):
                    cleaned_lines.append('')
        else:
            # Regular content - add as-is but handle spacing
            if line.strip() or (cleaned_lines and cleaned_lines[-1].strip()):
                cleaned_lines.append(line)
            i += 1
    
    return '\n'.join(cleaned_lines)
  
  
def _process_kable_tables(content: str) -> str:
    """
    Process R kable tables wrapped in cell-output-display blocks.
    
    Args:
        content: Markdown content containing kable table output
        
    Returns:
        str: Content with kable tables properly formatted
    """
    
    def replace_table_section(match):
        table_content = match.group(1)
        
        # If it's already a markdown table, keep as-is
        if '|' in table_content and '---' in table_content:
            return table_content.strip()
        
        # If it's a text table format, wrap in code block
        if '---' in table_content and ('  ' in table_content or '\n' in table_content):
            return f"```\n{table_content.strip()}\n```"
        
        # Default: wrap in code block
        return f"```\n{table_content.strip()}\n```"
    
    # Match ::: cell-output-display ... :::
    pattern = r'::: cell-output-display\s*\n(.*?)\n:::'
    content = re.sub(pattern, replace_table_section, content, flags=re.DOTALL)
    
    return content


def _process_table_captions(content: str) -> str:
    """
    Fix orphaned table captions and ensure proper code block closure.
    
    Args:
        content: Markdown content that may have caption formatting issues
        
    Returns:
        str: Content with properly formatted table captions
    """
    
    # Fix orphaned captions after code blocks
    def fix_table_caption(match):
        caption_text = match.group(1)
        return f"\n```\n\n**Table: {caption_text}**\n\n##"
    
    # Pattern: caption after content without proper code block closure
    pattern = r'\n\s*:\s+([^\n]+)\s*\n##'
    content = re.sub(pattern, fix_table_caption, content)
    
    # Fix incomplete table blocks with captions
    def fix_incomplete_table_block(match):
        table_content = match.group(1)
        caption_text = match.group(2)
        return f"{table_content}\n```\n\n**Table: {caption_text}**\n\n"
    
    # Pattern: table content + caption without proper closure
    pattern = r'(\n[^\n]*(?:---+[^\n]*\n[^\n]*)+)\s*:\s+([^\n]+)\s*\n'
    content = re.sub(pattern, fix_incomplete_table_block, content, flags=re.MULTILINE)
    
    return content

  
def _process_jupyter_display_tables(content: str) -> str:
    """
    Wrap unwrapped Jupyter text-based table displays in code blocks.
    
    Args:
        content: Markdown content containing Jupyter display output
        
    Returns:
        str: Content with text tables wrapped in code blocks
    """
    
    def wrap_text_table(match):
        table_text = match.group(1)
        return f"```\n{table_text.strip()}\n```"
    
    # Pattern for text tables:
    # - Header line 
    # - Separator line with dashes
    # - Data rows (but stop at empty lines or ```python/```r)
    # - Not already in code blocks
    pattern = r'(?<!```\n)(\s*\S+.*\n\s*[-=\s]{10,}\n(?:\s*\S+.*\n)*?)(?=\n\s*```|\n\s*$|\Z)'
    content = re.sub(pattern, wrap_text_table, content, flags=re.MULTILINE)
    
    return content
  
  
def _align_text_tables(content: str) -> str:
    """
    Fix alignment in text-based tables by properly spacing column headers.
    """
    logger.debug("align_text_tables called")
    
    lines = content.split('\n')
    code_blocks = _find_code_block_pairs(lines)
    
    if not code_blocks:
        logger.debug("No code block pairs found")
        return content
    
    logger.debug(f"Found {len(code_blocks)} code block pairs")
    
    for start_line, end_line in code_blocks:
        opening_line = lines[start_line].strip()
        
        # Skip python/r code blocks
        if opening_line in ['```python', '```r']:
            continue
        
        # Extract block content
        block_lines = lines[start_line + 1:end_line]
        if not block_lines:
            continue
        
        # Look for separator line with dashes
        separator_idx = -1
        for j, line in enumerate(block_lines):
            if line.count('-') >= 5:
                separator_idx = j
                break
        
        if separator_idx <= 0:  # No separator or separator is first line
            continue
        
        # Get header and separator lines
        header_line = block_lines[separator_idx - 1]
        separator_line = block_lines[separator_idx]
        
        # Parse header and find positions
        header_parts = re.split(r'\s{2,}|\t', header_line.strip())
        header_parts = [part.strip() for part in header_parts if part.strip()]
        
        dash_positions = _find_dash_positions(separator_line)
        
        # Determine column positions based on row names
        if len(header_parts) == len(dash_positions):
            positions = dash_positions  # No row names
        elif len(header_parts) == len(dash_positions) - 1:
            positions = dash_positions[1:]  # Has row names, skip first
        else:
            logger.warning(f"Cannot align table: {len(header_parts)} headers vs {len(dash_positions)} columns")
            continue
        
        # Align and replace header
        aligned_header = _align_header_to_positions(header_parts, positions, separator_line)
        lines[start_line + 1 + separator_idx - 1] = aligned_header
        
        logger.debug(f"Aligned table header: {header_parts} -> positions {positions}")
    
    return '\n'.join(lines)
  
  
def _remove_table_stylers(content: str) -> str:
    """
    Remove CSS style blocks from table output.
    
    Args:
        content: Markdown content containing CSS style blocks
        
    Returns:
        str: Content with CSS style blocks removed
    """
    
    # Remove <style> blocks but handle the case where ``` immediately follows
    # Pattern: <style>...</style>``` -> ```
    content = re.sub(r'<style[^>]*>.*?</style>\s*```', '```', content, flags=re.DOTALL | re.IGNORECASE)
    
    # Remove any remaining standalone style blocks
    content = re.sub(r'<style[^>]*>.*?</style>', '', content, flags=re.DOTALL | re.IGNORECASE)
    
    return content
  
  
def _find_code_block_pairs(lines: List[str]) -> List[Tuple[int, int]]:
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


def _find_dash_positions(separator_line: str) -> List[int]:
    """
    Find the starting positions of dash groups in a separator line.
    
    Returns:
        List of positions where dash groups start
    """
    positions = []
    in_dash_group = False
    
    for pos, char in enumerate(separator_line):
        if char == '-' and not in_dash_group:
            positions.append(pos)
            in_dash_group = True
        elif char != '-':
            in_dash_group = False
    
    return positions


def _align_header_to_positions(header_parts: List[str], positions: List[int], separator_line: str) -> str:
    """
    Align header parts to their corresponding dash positions.
    
    Returns:
        Aligned header string
    """
    new_header = [' '] * len(separator_line)
    
    for header, pos in zip(header_parts, positions):
        for i, char in enumerate(header):
            if pos + i < len(new_header):
                new_header[pos + i] = char
    
    return ''.join(new_header).rstrip()
