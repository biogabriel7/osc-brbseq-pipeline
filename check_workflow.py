#!/usr/bin/env python3
import os
import sys
import subprocess

def main():
    """Check the Snakemake workflow for errors."""
    print("Checking Snakemake workflow...")
    
    # Run snakemake with dry-run and print shell commands
    cmd = ["snakemake", "-n", "-p"]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        # Check if there were any errors
        if result.returncode != 0:
            print("ERROR: Snakemake workflow has errors:")
            print(result.stderr)
            return 1
        
        # Print the workflow plan
        print("Workflow plan:")
        print(result.stdout)
        
        print("\nNo errors found in the workflow!")
        return 0
        
    except Exception as e:
        print(f"ERROR: Failed to run Snakemake: {e}")
        return 1

if __name__ == "__main__":
    sys.exit(main()) 