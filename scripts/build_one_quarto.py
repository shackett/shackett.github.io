#!/usr/bin/env python3
import sys
import os
import subprocess
from pathlib import Path
import tempfile
import shutil
import json

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

# Rest of the code stays the same, but replace the yaml.dump line with:
def build_one_qmd(input_file, output_file):
    """Render a single .qmd file to Jekyll-compatible markdown"""
    
    print(f"* rendering {input_file}")
    
    # Get the base name for figure paths (matching your R setup)
    stem = Path(input_file).stem
    base_name = stem.replace("_source/", "").replace("_source\\", "")
    
    # Load base config and customize for this file
    config = load_quarto_config()
    config['format']['gfm']['fig-path'] = f"figure/{base_name}/"
    config['execute']['cache-path'] = f"cache/{base_name}/"
    
    # Create a temporary directory for this render
    temp_dir = tempfile.mkdtemp()
    temp_yml = os.path.join(temp_dir, "_quarto.yml")
    
    # Write the customized config as YAML
    with open(temp_yml, 'w') as f:
        f.write(config_to_yaml_string(config))
    
    # Create a temporary copy of the input file in temp_dir
    temp_input = os.path.join(temp_dir, os.path.basename(input_file))
    shutil.copy2(input_file, temp_input)
    
    try:
        # Render with Quarto (let it output in temp_dir, then move)
        cmd = [
            "quarto", "render", temp_input,
            "--to", "gfm"
        ]
        
        result = subprocess.run(cmd, cwd=temp_dir, check=True, capture_output=True, text=True)
        
        # Find the output file (should be .md in temp_dir)
        temp_output = os.path.join(temp_dir, Path(input_file).stem + ".md")
        
        if os.path.exists(temp_output):
            # Move to final destination
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
            shutil.move(temp_output, output_file)
            clean_jekyll_output(output_file)
        else:
            raise FileNotFoundError(f"Expected output file not found: {temp_output}")
            
    except subprocess.CalledProcessError as e:
        if os.path.exists(output_file):
            os.unlink(output_file)
        print(f"Error: {e}")
        print(f"Stdout: {e.stdout}")
        print(f"Stderr: {e.stderr}")
        sys.exit(1)
    finally:
        # Clean up temp directory
        shutil.rmtree(temp_dir, ignore_errors=True)

def clean_jekyll_output(file_path):
    """Clean up the output file for Jekyll compatibility"""
    with open(file_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Remove any Quarto-specific artifacts that might interfere with Jekyll
    # You can add more cleaning rules here as needed
    
    with open(file_path, 'w', encoding='utf-8') as f:
        f.write(content)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python build_one_quarto.py input.qmd output.md")
        sys.exit(1)
    
    build_one_qmd(sys.argv[1], sys.argv[2])
