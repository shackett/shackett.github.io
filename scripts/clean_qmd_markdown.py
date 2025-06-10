"""
Quarto Markdown Cleaning Module

This module contains functions to clean and process Quarto-generated markdown
to make it compatible with Jekyll.
"""

import re
import logging
from pathlib import Path


def remove_quarto_artifacts(content: str) -> str:
    """
    Remove basic Quarto cell artifacts.
    
    Args:
        content: Raw markdown content from Quarto
        
    Returns:
        str: Content with Quarto artifacts removed
    """
    """Remove basic Quarto cell artifacts"""
    
    # Replace cell dividers and convert to Jekyll format
    content = re.sub(r'::+\s*cell\s*\n', '', content)
    content = re.sub(r'::+\s*\n', '', content)
    
    # Convert ``` {.python .cell-code} to ```python
    content = re.sub(r'``` \{\.python \.cell-code\}', '```python', content)
    content = re.sub(r'``` \{\.r \.cell-code\}', '```r', content)
    
    # Remove cell output wrappers (including stderr)
    content = re.sub(r'``` \{\.cell-output \.cell-output-stdout\}\s*\n', '', content)
    content = re.sub(r'``` \{\.cell-output \.cell-output-stderr\}\s*\n', '', content)
    content = re.sub(r'::: \{\.cell-output \.cell-output-stdout\}\s*\n', '', content)
    content = re.sub(r'::: \{\.cell-output \.cell-output-stderr\}\s*\n', '', content)
    
    return content


def process_quarto_table(content: str) -> str:
    """
    Process Quarto table output for Jekyll compatibility.
    
    Handles kable tables wrapped in cell-output-display blocks and
    fixes orphaned table captions.
    
    Args:
        content: Markdown content containing Quarto table output
        
    Returns:
        str: Content with properly formatted tables
    """
    
    # Handle kable tables wrapped in cell-output-display
    def replace_table_section(match):
        table_content = match.group(1)
        
        # If it's a markdown table (contains | and ---), keep as-is
        if '|' in table_content and '---' in table_content:
            return table_content.strip()
        
        # Check if it's a text table format (with dashes for table structure)
        if '---' in table_content and ('  ' in table_content or '\n' in table_content):
            # It's likely a text-formatted table, wrap in code block
            return f"```\n{table_content.strip()}\n```"
        
        # Otherwise, wrap in code block for now
        return f"```\n{table_content.strip()}\n```"
    
    # Pattern to match ::: cell-output-display ... :::
    pattern = r'::: cell-output-display\s*\n(.*?)\n:::'
    content = re.sub(pattern, replace_table_section, content, flags=re.DOTALL)
    
    # Fix orphaned table captions that appear after code blocks
    # The key issue: ensure code block is closed AND caption is formatted properly
    def fix_table_caption(match):
        caption_text = match.group(1)
        # Ensure code block is properly closed on its own line, add caption, then continue
        return f"\n```\n\n**Table: {caption_text}**\n\n##"
    
    # Look for pattern where table content ends but no closing ``` before caption
    # This handles the specific case you're seeing
    pattern = r'\n\s*:\s+([^\n]+)\s*\n##'
    content = re.sub(pattern, fix_table_caption, content)
    
    # Also handle cases where we have incomplete code blocks before captions
    # Look for table content followed by caption (without proper closing)
    def fix_incomplete_table_block(match):
        table_content = match.group(1)
        caption_text = match.group(2)
        
        # Ensure the table block is properly closed on its own line and caption is formatted
        return f"{table_content}\n```\n\n**Table: {caption_text}**\n\n"
    
    # Pattern for table content + caption without proper code block closure
    pattern = r'(\n[^\n]*(?:---+[^\n]*\n[^\n]*)+)\s*:\s+([^\n]+)\s*\n'
    content = re.sub(pattern, fix_incomplete_table_block, content, flags=re.MULTILINE)
    
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
    logger = logging.getLogger(__name__)
    
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
