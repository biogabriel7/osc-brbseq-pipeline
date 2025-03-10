#!/usr/bin/env python3
"""
Script to fix syntax errors in Snakefile by adding semicolons after rule declarations.
This addresses the common "Statements must be separated by newlines or semicolons" errors.
"""

import re
import sys
import os
import shutil
from datetime import datetime

def fix_snakefile(file_path):
    """
    Fix the Snakefile syntax by adding semicolons after rule declarations.
    
    Args:
        file_path: Path to the Snakefile to fix
    
    Returns:
        True if changes were made, False otherwise
    """
    # Create a backup of the original file
    backup_path = f"{file_path}.bak.{datetime.now().strftime('%Y%m%d%H%M%S')}"
    shutil.copy2(file_path, backup_path)
    print(f"Created backup at {backup_path}")
    
    # Read the file content
    with open(file_path, 'r') as f:
        content = f.read()
    
    # Fix rule declarations by adding semicolons
    # Pattern: matches "rule name:" not followed by a semicolon or newline+whitespace+input/output/etc.
    pattern = r'(rule\s+[a-zA-Z0-9_]+:)(?!\s*;)(?!\s*\n\s+(input|output|log|params|threads|resources|shell|run|script|conda|container|message|priority|benchmark|version|envmodules|shadow|group|cache|handover|default_target|localrule))'
    
    # Also fix ruleorder declarations
    ruleorder_pattern = r'(ruleorder:.*?)(?!\s*;)(?=\s*\n)'
    
    # Apply the fixes
    fixed_content = re.sub(pattern, r'\1;', content)
    fixed_content = re.sub(ruleorder_pattern, r'\1;', fixed_content)
    
    # Fix onsuccess and onerror declarations
    onsuccess_pattern = r'(onsuccess:)(?!\s*;)(?=\s*\n)'
    onerror_pattern = r'(onerror:)(?!\s*;)(?=\s*\n)'
    
    fixed_content = re.sub(onsuccess_pattern, r'\1;', fixed_content)
    fixed_content = re.sub(onerror_pattern, r'\1;', fixed_content)
    
    # Check if any changes were made
    if fixed_content == content:
        print("No changes needed in the file.")
        return False
    
    # Write the fixed content back to the file
    with open(file_path, 'w') as f:
        f.write(fixed_content)
    
    print(f"Fixed syntax errors in {file_path}")
    return True

if __name__ == "__main__":
    if len(sys.argv) > 1:
        snakefile_path = sys.argv[1]
    else:
        snakefile_path = "Snakefile"  # Default path
    
    if not os.path.exists(snakefile_path):
        print(f"Error: File {snakefile_path} not found.")
        sys.exit(1)
    
    success = fix_snakefile(snakefile_path)
    if success:
        print("Syntax errors fixed. Please check the file and rerun your workflow.")
    else:
        print("No syntax errors found or fixed.") 