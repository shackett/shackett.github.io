#!/usr/bin/env python3
"""
Quarto Document Builder for Jekyll

This module renders Quarto (.qmd) documents to Jekyll-compatible markdown.
It handles both R and Python code execution via knitr engine.
"""

import sys
import os
import subprocess
import tempfile
import shutil
import json
import logging
from pathlib import Path
from typing import Dict, Any

from utils import setup_logging
from clean_qmd_markdown import clean_jekyll_output


def load_quarto_config() -> Dict[str, Any]:
    """
    Load the base Quarto configuration from JSON file.
    
    Returns:
        Dict[str, Any]: Quarto configuration dictionary
        
    Raises:
        FileNotFoundError: If quarto_config.json is not found
        json.JSONDecodeError: If JSON is malformed
    """
    config_path = Path(__file__).parent / "quarto_config.json"
    with open(config_path, 'r', encoding='utf-8') as f:
        return json.load(f)


def config_to_yaml_string(config: Dict[str, Any]) -> str:
    """
    Convert config dictionary to YAML string manually.
    
    This avoids the pyyaml dependency by implementing basic YAML serialization.
    
    Args:
        config: Configuration dictionary to convert
        
    Returns:
        str: YAML-formatted string
    """
    def dict_to_yaml(d: Dict[str, Any], indent: int = 0) -> str:
        """Recursively convert dictionary to YAML format."""
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


def copy_figures_to_jekyll_location(temp_dir: str, filename_stem: str) -> None:
    """
    Copy generated figures from temporary directory to Jekyll location.
    
    Quarto generates figures in various locations within the temp directory.
    This function searches for all figure files and copies them to the
    appropriate Jekyll figure directory.
    
    Args:
        temp_dir: Path to temporary rendering directory
        filename_stem: Base filename without extension (used for figure path)
    """
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
    
    # Check primary figure location
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
    
    # Check alternative files location
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
    
    # Search entire temp directory for any other figure files
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


def build_one_qmd(input_file: str, output_file: str) -> None:
    """
    Render a single .qmd file to Jekyll-compatible markdown.
    
    This function:
    1. Creates a temporary directory for rendering
    2. Configures Quarto with appropriate figure and cache paths
    3. Renders the document using Quarto
    4. Copies generated figures to Jekyll location
    5. Cleans the output markdown for Jekyll compatibility
    
    Args:
        input_file: Path to input .qmd file
        output_file: Path where output .md file should be saved
        
    Raises:
        subprocess.CalledProcessError: If Quarto rendering fails
        FileNotFoundError: If expected output file is not generated
    """
    logger = logging.getLogger(__name__)
    
    logger.info(f"Rendering {input_file}")
    
    # Get the base name for figure paths (matching R setup)
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
    
    with open(temp_yml, 'w', encoding='utf-8') as f:
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
        subprocess.run(cmd, cwd=temp_dir, check=True, capture_output=True, text=True)
        logger.debug(f"Files in temp_dir after render: {os.listdir(temp_dir)}")
        
        if os.path.exists(temp_output):
            # Debug: show what Quarto produced
            with open(temp_output, 'r', encoding='utf-8') as f:
                raw_content = f.read()
                logger.debug(f"Raw Quarto output (first 1000 chars):\n{raw_content[:1000]}")
            
            # Copy figures to Jekyll location BEFORE moving files
            logger.info("Copying figures before moving output file...")
            copy_figures_to_jekyll_location(temp_dir, filename_stem)
            
            # Move to final destination
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
            shutil.move(temp_output, output_file)
            
            # Clean the output file using the cleaning module
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
        # Clean up temp directory
        shutil.rmtree(temp_dir, ignore_errors=True)


def main() -> None:
    """Main entry point for the script."""
    setup_logging()
    
    if len(sys.argv) != 3:
        print("Usage: python build_one_quarto.py input.qmd output.md")
        sys.exit(1)
    
    build_one_qmd(sys.argv[1], sys.argv[2])


if __name__ == "__main__":
    main()
