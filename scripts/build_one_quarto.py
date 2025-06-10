#!/usr/bin/env python3
import sys
import os
import subprocess
from pathlib import Path
import tempfile
import shutil
import json
import logging
from utils import setup_logging

def load_quarto_config():
    """Load the base Quarto configuration from JSON file"""
    config_path = Path(__file__).parent / "quarto_config.json"
    with open(config_path, 'r') as f:
        return json.load(f)

def config_to_yaml_string(config):
    """Convert config dict to YAML string manually (avoiding pyyaml dependency)"""
    def dict_to_yaml(d, indent=0):
        yaml_str = ""
        for key, value in d.items():
            spaces = "  " * indent
            if isinstance(value, dict):
                yaml_str += f"{spaces}{key}:\n"
                yaml_str += dict_to_yaml(value, indent + 1)
            elif isinstance(value, bool):
                yaml_str += f"{spaces}{key}: {str(value).lower()}\n"
            else:
                yaml_str += f"{spaces}{key}: {value}\n"
        return yaml_str
    
    return dict_to_yaml(config)

def process_quarto_table(content):
    """Process Quarto table output for Jekyll compatibility"""
    import re
    
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
    
    # Fix orphaned table captions (lines starting with ": " that are alone)
    content = re.sub(r'\n\s*:\s+([^\n]+)\n##', r'\n\n**\1**\n\n##', content)
    
    return content

def process_quarto_figures(content, filename_stem):
    """Process figure paths and display blocks"""
    import re
    
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
    
    # Convert absolute figure paths to Jekyll-relative paths
    # The key insight: Jekyll serves from the root, so /figure/... is correct
    # But our markdown has figure/... (relative), which Jekyll interprets relative to the post URL
    # We need to change figure/source/... to /figure/source/...
    content = re.sub(r'!\[([^\]]*)\]\(figure/source/', r'![\1](/figure/source/', content)
    
    return content

def process_code_blocks_and_output(content):
    """Process code blocks and their output with proper spacing"""
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

def remove_quarto_artifacts(content):
    """Remove basic Quarto cell artifacts"""
    import re
    
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

def copy_figures_to_jekyll_location(temp_dir, filename_stem):
    """Copy generated figures from temp directory to Jekyll location"""
    import shutil
    import os
    
    logger = logging.getLogger(__name__)
    
    # Source paths in temp directory
    temp_figure_dir = os.path.join(temp_dir, "figure", "source", filename_stem)
    temp_files_dir = os.path.join(temp_dir, f"{filename_stem}_files", "figure-markdown")
    
    # Destination path in Jekyll
    dest_figure_dir = f"figure/source/{filename_stem}"
    
    # Create destination directory
    os.makedirs(dest_figure_dir, exist_ok=True)
    
    # Copy from both possible locations
    figures_copied = 0
    
    logger.debug(f"Checking temp figure directory: {os.path.abspath(temp_figure_dir)}")
    logger.debug(f"Checking temp files directory: {os.path.abspath(temp_files_dir)}")
    logger.debug(f"Destination directory: {os.path.abspath(dest_figure_dir)}")
    
    if os.path.exists(temp_figure_dir):
        logger.debug(f"Found temp figure dir: {temp_figure_dir}")
        files_in_dir = os.listdir(temp_figure_dir)
        logger.debug(f"Files in temp figure dir: {files_in_dir}")
        
        for file in files_in_dir:
            if file.endswith(('.png', '.jpg', '.jpeg', '.svg', '.pdf')):
                src = os.path.join(temp_figure_dir, file)
                dst = os.path.join(dest_figure_dir, file)
                shutil.copy2(src, dst)
                figures_copied += 1
                logger.info(f"Copied figure: {src} -> {dst}")
    else:
        logger.debug(f"Temp figure directory does not exist: {temp_figure_dir}")
    
    if os.path.exists(temp_files_dir):
        logger.debug(f"Found temp files dir: {temp_files_dir}")
        files_in_dir = os.listdir(temp_files_dir)
        logger.debug(f"Files in temp files dir: {files_in_dir}")
        
        for file in files_in_dir:
            if file.endswith(('.png', '.jpg', '.jpeg', '.svg', '.pdf')):
                src = os.path.join(temp_files_dir, file)
                dst = os.path.join(dest_figure_dir, file)
                shutil.copy2(src, dst)
                figures_copied += 1
                logger.info(f"Copied figure: {src} -> {dst}")
    else:
        logger.debug(f"Temp files directory does not exist: {temp_files_dir}")
    
    # Also check for any other figure files in the temp directory
    logger.debug("Searching for all figure files in temp directory...")
    for root, dirs, files in os.walk(temp_dir):
        for file in files:
            if file.endswith(('.png', '.jpg', '.jpeg', '.svg', '.pdf')):
                full_path = os.path.join(root, file)
                logger.debug(f"Found figure file: {full_path}")
                
                # Copy to destination if not already copied
                dst = os.path.join(dest_figure_dir, file)
                if not os.path.exists(dst):
                    shutil.copy2(full_path, dst)
                    figures_copied += 1
                    logger.info(f"Copied additional figure: {full_path} -> {dst}")
    
    logger.info(f"Total figures copied: {figures_copied} to {os.path.abspath(dest_figure_dir)}")
    
    # List what's actually in the destination
    if os.path.exists(dest_figure_dir):
        dest_files = os.listdir(dest_figure_dir)
        logger.info(f"Files now in destination: {dest_files}")
    else:
        logger.warning(f"Destination directory was not created: {dest_figure_dir}")

def clean_jekyll_output(file_path):
    """Clean up the output file for Jekyll compatibility"""
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

def build_one_qmd(input_file, output_file):
    """Render a single .qmd file to Jekyll-compatible markdown"""
    logger = logging.getLogger(__name__)
    
    logger.info(f"Rendering {input_file}")
    
    # Get the base name for figure paths (matching your R setup)
    filename_stem = Path(input_file).stem
    base_name = f"source/{filename_stem}"
    
    # Load base config and customize for this file
    config = load_quarto_config()
    config['format']['markdown']['fig-path'] = f"figure/{base_name}/"
    config['execute']['cache-path'] = f"cache/{base_name}/"
    
    # Create a temporary directory for this render
    temp_dir = tempfile.mkdtemp()
    temp_yml = os.path.join(temp_dir, "_quarto.yml")
    temp_output = os.path.join(temp_dir, Path(input_file).stem + ".md")
    
    logger.debug(f"Using temp directory: {temp_dir}")
    
    # Write the customized config as YAML
    yaml_config = config_to_yaml_string(config)
    logger.debug(f"Generated YAML config:\n{yaml_config}")
    
    with open(temp_yml, 'w') as f:
        f.write(yaml_config)
    
    # Create a temporary copy of the input file in temp_dir
    temp_input = os.path.join(temp_dir, os.path.basename(input_file))
    shutil.copy2(input_file, temp_input)
    
    try:
        # Use markdown format and ensure execution is enabled
        cmd = [
            "quarto", "render", temp_input,
            "--to", "markdown",
            "--execute"
        ]
        
        logger.debug(f"Running: {' '.join(cmd)}")
        result = subprocess.run(cmd, cwd=temp_dir, check=True, capture_output=True, text=True)
        
        logger.debug(f"Quarto stdout: {result.stdout}")
        logger.debug(f"Files in temp_dir after render: {os.listdir(temp_dir)}")
        
        if os.path.exists(temp_output):
            # Debug: show what Quarto produced
            with open(temp_output, 'r') as f:
                raw_content = f.read()
                logger.debug(f"Raw Quarto output (first 1000 chars):\n{raw_content[:1000]}")
            
            # Copy figures to Jekyll location BEFORE moving files
            logger.info("Copying figures before moving output file...")
            copy_figures_to_jekyll_location(temp_dir, filename_stem)
            
            # Move to final destination
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
            shutil.move(temp_output, output_file)
            
            # Clean the output file
            clean_jekyll_output(output_file)
            
            logger.debug(f"Successfully built {input_file}")
        else:
            raise FileNotFoundError(f"Expected output file not found: {temp_output}")
            
    except subprocess.CalledProcessError as e:
        logger.error(f"Quarto render failed for {input_file}")
        logger.error(f"Command: {' '.join(cmd)}")
        logger.error(f"Return code: {e.returncode}")
        if e.stdout:
            logger.error(f"Stdout: {e.stdout}")
        if e.stderr:
            logger.error(f"Stderr: {e.stderr}")
        
        # Clean up failed output file
        if os.path.exists(output_file):
            os.unlink(output_file)
        sys.exit(1)
    finally:
        # Temporarily disable cleanup for debugging
        logger.info(f"DEBUGGING: Temp directory preserved at: {temp_dir}")
        # shutil.rmtree(temp_dir, ignore_errors=True)

if __name__ == "__main__":
    setup_logging()
    
    if len(sys.argv) != 3:
        print("Usage: python build_one_quarto.py input.qmd output.md")
        sys.exit(1)
    
    build_one_qmd(sys.argv[1], sys.argv[2])
