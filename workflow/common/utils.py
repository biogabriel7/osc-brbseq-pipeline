#!/usr/bin/env python3
"""
Utility functions for RNA-seq Snakemake workflow
"""

import os
import csv

def get_samples_from_metasheet(config_file, test_mode=False):
    """
    Parse metasheet and extract sample information.
    
    Args:
        config_file: Dictionary containing workflow configuration
        test_mode: If True, only return the first sample for testing
        
    Returns:
        Dictionary of samples with read paths and metadata
    """
    samples = {}
    try:
        with open(config_file["metasheet_path"], "r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                sample_name = row["sample"]
                samples[sample_name] = {
                    "R1": row["R1"],
                    "R2": row["R2"],
                    "srr": os.path.basename(row["R1"]).replace("_R1.fastq.gz", "")  # Get SRR ID from filename
                }
                print(f"Loaded sample: {sample_name}, SRR: {samples[sample_name]['srr']}")
                
                # If in test mode, just take the first sample and break
                if test_mode and len(samples) == 1:
                    print("TEST MODE: Using only the first sample")
                    break
                    
    except FileNotFoundError:
        print(f"Error: Metasheet not found at {config_file['metasheet_path']}")
        raise
    except KeyError as e:
        print(f"Error: Missing required column in metasheet: {e}")
        raise
        
    return samples

def create_workflow_dirs(base_dirs=None):
    """
    Create all necessary directories for the workflow.
    
    Args:
        base_dirs: Optional list of additional directories to create
    """
    if base_dirs is None:
        base_dirs = []
        
    # Basic workflow directory structure
    dirs = [
        # Log directories
        "logs/fastqc", "logs/multiqc", "logs/trimming", "logs/star", 
        "logs/workflow", "logs/slurm",
        
        # Analysis directories
        "Analysis/QC/FastQC", "Analysis/QC/MultiQC", "Analysis/QC/Trimming",
        "Analysis/QC/Trimming/Reports", "Analysis/QC/Trimming/MultiQC", 
        "Analysis/Trimmed", "Analysis/Alignment/STAR", "Analysis/Alignment/MultiQC",
        
        # Resource directories
        "resources/metadata"
    ]
    
    # Add any additional directories
    dirs.extend(base_dirs)
    
    # Create all directories
    for dir_path in dirs:
        os.makedirs(dir_path, exist_ok=True)
        print(f"Created directory: {dir_path}")