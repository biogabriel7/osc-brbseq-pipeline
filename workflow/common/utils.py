#!/usr/bin/env python3
"""
Utility functions for RNA-seq Snakemake workflow
"""

import os
import csv
import pandas as pd

def get_samples_from_metasheet(config, test_mode=False):
    """
    Read sample information from metasheet
    
    Args:
        config (dict): Configuration dictionary
        test_mode (bool): If True, only return the first sample
        
    Returns:
        dict: Dictionary of sample information
    """
    metasheet_path = config.get("metasheet_path", "resources/metadata/samples.csv")
    
    if not os.path.exists(metasheet_path):
        raise FileNotFoundError(f"Metasheet file not found: {metasheet_path}")
    
    # Read metasheet
    try:
        df = pd.read_csv(metasheet_path)
        required_columns = ["sample", "R1", "R2"]
        for col in required_columns:
            if col not in df.columns:
                raise ValueError(f"Required column '{col}' not found in metasheet")
        
        # Convert to dictionary
        samples = {}
        for _, row in df.iterrows():
            sample_name = row["sample"]
            samples[sample_name] = {
                "R1": row["R1"],
                "R2": row["R2"]
            }
            
        # If in test mode, only return the first sample
        if test_mode and samples:
            first_sample = list(samples.keys())[0]
            return {first_sample: samples[first_sample]}
            
        return samples
    except Exception as e:
        raise ValueError(f"Error reading metasheet: {str(e)}")

def create_workflow_dirs():
    """
    Create all necessary directories for the workflow
    """
    # Create main analysis directories
    os.makedirs("Analysis/QC/FastQC", exist_ok=True)
    os.makedirs("Analysis/QC/MultiQC", exist_ok=True)
    os.makedirs("Analysis/QC/Trimming", exist_ok=True)
    os.makedirs("Analysis/QC/Trimming/Reports", exist_ok=True)
    os.makedirs("Analysis/QC/Trimming/MultiQC", exist_ok=True)
    os.makedirs("Analysis/Trimmed", exist_ok=True)
    os.makedirs("Analysis/Alignment/STAR", exist_ok=True)
    os.makedirs("Analysis/Alignment/MultiQC", exist_ok=True)
    os.makedirs("Analysis/Counts/FeatureCounts", exist_ok=True)
    os.makedirs("Analysis/Counts/MultiQC", exist_ok=True)
    
    # Create log directories
    os.makedirs("logs/fastqc", exist_ok=True)
    os.makedirs("logs/multiqc", exist_ok=True)
    os.makedirs("logs/trimming", exist_ok=True)
    os.makedirs("logs/star", exist_ok=True)
    os.makedirs("logs/counts", exist_ok=True)
    os.makedirs("logs/workflow", exist_ok=True)
    
    # Create resource directories
    os.makedirs("resources/metadata", exist_ok=True)
    os.makedirs("resources/genome", exist_ok=True)
    os.makedirs("resources/config", exist_ok=True)