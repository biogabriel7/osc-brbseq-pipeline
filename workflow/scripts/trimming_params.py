#!/usr/bin/env python3
"""
Script to analyze MultiQC FastQC results and generate sample-specific trimming parameters.
This script reads multiqc_fastqc.txt, analyzes quality metrics for each sample,
and generates a JSON file with sample-specific trimming parameters.
Modified to work with fastp rather than Trimmomatic.

Usage:
    python trimming_params.py --input multiqc_fastqc.txt --output trimming_params.json [--config params.yaml]

Author: Gabriel D.
Date: 2025-03-05
"""

import argparse
import csv
import json
import yaml
import os
from collections import defaultdict

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Generate trimming parameters from MultiQC FastQC results')
    parser.add_argument('--input', required=True, help='Input multiqc_fastqc.txt file')
    parser.add_argument('--output', required=True, help='Output JSON file for trimming parameters')
    parser.add_argument('--config', help='Snakemake config YAML file with default parameters')
    return parser.parse_args()

def load_default_params(config_file):
    """Load default trimming parameters from config file."""
    if not config_file or not os.path.exists(config_file):
        # Default parameters if config file not provided
        return {
            "leading_quality": 3,
            "trailing_quality": 3,
            "sliding_window": "4:15",
            "min_length": 36
        }
    
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    
    # Extract trimming parameters from config
    params = {
        "leading_quality": config.get("leading_quality", 3),
        "trailing_quality": config.get("trailing_quality", 3),
        "sliding_window": config.get("sliding_window", "4:15"),
        "min_length": config.get("min_length", 36)
    }
    
    return params

# This is a partial update for the read_fastqc_data function in trimming_params.py

def read_fastqc_data(input_file):
    """Read and parse the multiqc_fastqc.txt file."""
    samples_data = {}
    
    try:
        with open(input_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                # Check if 'Sample' column exists, otherwise try 'Filename'
                if 'Sample' in row:
                    sample_name = row['Sample']
                elif 'Filename' in row:
                    # Extract sample name from filename
                    filename = row['Filename']
                    # Remove _fastqc.zip suffix if present
                    sample_name = filename.replace('_fastqc.zip', '')
                    # Remove .fastq.gz suffix if present
                    sample_name = sample_name.replace('.fastq.gz', '')
                else:
                    # Skip row if no identifiable sample name
                    continue
                
                # Extract read direction (R1 or R2)
                read_direction = 'R1' if '_R1' in sample_name else 'R2' if '_R2' in sample_name else 'unknown'
                
                # Extract SRR ID for more reliable matching with Snakefile
                # This should identify the SRR ID pattern in the filename
                import re
                srr_match = re.search(r'(SRR\d+)', sample_name)
                if srr_match:
                    base_sample = srr_match.group(1)
                else:
                    # Extract base sample name without read direction
                    base_sample = sample_name
                    if '_R1' in sample_name:
                        base_sample = sample_name.replace('_R1', '')
                    elif '_R2' in sample_name:
                        base_sample = sample_name.replace('_R2', '')
                
                # Store data by base sample and read direction
                if base_sample not in samples_data:
                    samples_data[base_sample] = {}
                
                # Store the row with original case preserved
                samples_data[base_sample][read_direction] = row
        
        # Debug output
        print(f"Processed {len(samples_data)} samples from {input_file}")
        for sample in samples_data:
            print(f"Found sample: {sample}")
        
        if not samples_data:
            print("Warning: No samples were found in the input file")
            # Print first few lines of the file for debugging
            with open(input_file, 'r') as f:
                print("File content preview:")
                for i, line in enumerate(f):
                    if i < 5:  # Print first 5 lines
                        print(line.strip())
                    else:
                        break
    
    except Exception as e:
        print(f"Error reading FastQC data: {str(e)}")
        raise
    
    return samples_data

def analyze_sample(sample_name, sample_data, default_params):
    """
    Analyze FastQC metrics for a sample and determine appropriate trimming parameters.
    
    Args:
        sample_name: Name of the sample
        sample_data: Dict containing R1 and R2 FastQC metrics
        default_params: Dict containing default trimming parameters
    
    Returns:
        Dict with sample-specific trimming parameters
    """
    # Start with default parameters
    params = default_params.copy()
    
    # Determine if any special trimming is needed based on QC metrics
    r1_data = sample_data.get('R1', {})
    r2_data = sample_data.get('R2', {})
    
    # Flag to track if any customization was applied
    customized = False
    
    # Check for adapter content issues - case-insensitive search
    adapter_key = next((k for k in r1_data.keys() if k.lower() == 'adapter_content'), None)
    if adapter_key:
        if (r1_data.get(adapter_key) in ['warn', 'fail'] or 
            r2_data.get(adapter_key) in ['warn', 'fail']):
            # With fastp, we're using auto-detection, but we could adjust stringency if needed
            params['adapter_stringency'] = 0.8  # For reference - fastp handles this differently
            customized = True
    
    # Check for quality issues - case-insensitive search
    quality_key = next((k for k in r1_data.keys() if k.lower() == 'per_base_sequence_quality'), None)
    if quality_key:
        if (r1_data.get(quality_key) in ['warn', 'fail'] or 
            r2_data.get(quality_key) in ['warn', 'fail']):
            # More aggressive quality trimming
            params['leading_quality'] = max(5, params['leading_quality'])
            params['trailing_quality'] = max(5, params['trailing_quality'])
            params['sliding_window'] = "4:20"  # Higher quality threshold
            customized = True
    
    # Check for sequence length distribution issues
    length_key = next((k for k in r1_data.keys() if k.lower() == 'sequence_length_distribution'), None)
    if length_key:
        if (r1_data.get(length_key) in ['warn', 'fail'] or 
            r2_data.get(length_key) in ['warn', 'fail']):
            # Adjust minimum length requirement
            # Get average sequence length, ensure it's a number
            seq_length_key = next((k for k in r1_data.keys() if k.lower() == 'avg_sequence_length'), None)
            
            if seq_length_key and r1_data.get(seq_length_key):
                try:
                    # Try to convert to float first
                    seq_length = float(r1_data.get(seq_length_key, 150))
                    # Then convert to int for calculation
                    params['min_length'] = max(36, int(seq_length * 0.7))  # At least 70% of average length
                except (ValueError, TypeError):
                    # If conversion fails, use default value
                    print(f"Warning: Could not convert sequence length '{r1_data.get(seq_length_key)}' to number. Using default.")
                    params['min_length'] = 36
            else:
                # If no length info, use default
                params['min_length'] = 36
                
            customized = True
    
    # Check for N content
    n_content_key = next((k for k in r1_data.keys() if k.lower() == 'per_base_n_content'), None)
    if n_content_key:
        if (r1_data.get(n_content_key) in ['warn', 'fail'] or 
            r2_data.get(n_content_key) in ['warn', 'fail']):
            # Max N filter - fastp has n_base_limit parameter
            params['n_base_limit'] = 5  # Limit bases with N
            customized = True
    
    # Check for duplication levels
    duplication_key = next((k for k in r1_data.keys() if k.lower() == 'sequence_duplication_levels'), None)
    dedup_pct_key = next((k for k in r1_data.keys() if k.lower() == 'total_deduplicated_percentage'), None)
    
    if duplication_key and dedup_pct_key:
        if (r1_data.get(duplication_key) in ['warn', 'fail'] or 
            r2_data.get(duplication_key) in ['warn', 'fail']):
            # If it's very high, might suggest PCR duplicates
            try:
                total_dedup_pct = float(r1_data.get(dedup_pct_key, 100))
                if total_dedup_pct < 50:  # High duplication rate
                    params['deduplication'] = True
                    customized = True
            except (ValueError, TypeError):
                print(f"Warning: Could not convert deduplication percentage '{r1_data.get(dedup_pct_key)}' to number.")
    
    # Check for overrepresented sequences
    overrep_key = next((k for k in r1_data.keys() if k.lower() == 'overrepresented_sequences'), None)
    if overrep_key:
        if (r1_data.get(overrep_key) in ['warn', 'fail'] or 
            r2_data.get(overrep_key) in ['warn', 'fail']):
            # Use fastp's low complexity filter to help with this
            params['low_complexity_filter'] = True
            customized = True
    
    # Check for per base sequence content issues (often related to biases)
    base_content_key = next((k for k in r1_data.keys() if k.lower() == 'per_base_sequence_content'), None)
    if base_content_key:
        if (r1_data.get(base_content_key) in ['warn', 'fail'] or 
            r2_data.get(base_content_key) in ['warn', 'fail']):
            # This could be due to various factors
            # For NextSeq/NovaSeq, it could be poly-G tails
            # Let's check if the sample might be from these platforms by name pattern
            if "NextSeq" in sample_name or "NovaSeq" in sample_name:
                params['trim_poly_g'] = True
                customized = True
            # Alternatively, check if GC content is very skewed
            gc_key = next((k for k in r1_data.keys() if k.lower() == '%gc'), None)
            if gc_key:
                try:
                    gc_content = float(r1_data.get(gc_key, 50))
                    if gc_content > 60 or gc_content < 40:
                        params['correction'] = True  # Enable correction
                        customized = True
                except (ValueError, TypeError):
                    pass
    
    # Tag whether this sample needed customization
    params['customized'] = customized
    
    return params

def generate_trimming_params(samples_data, default_params):
    """Generate sample-specific trimming parameters based on QC metrics."""
    trimming_params = {}
    
    for sample_name, sample_data in samples_data.items():
        # Analyze sample and get appropriate parameters
        params = analyze_sample(sample_name, sample_data, default_params)
        
        # Store parameters for this sample
        trimming_params[sample_name] = params
    
    return trimming_params

def main():
    """Main function to run the script."""
    args = parse_args()
    
    # Load default parameters
    default_params = load_default_params(args.config)
    
    # Read FastQC data
    samples_data = read_fastqc_data(args.input)
    
    # Generate trimming parameters
    trimming_params = generate_trimming_params(samples_data, default_params)
    
    # Output JSON file
    with open(args.output, 'w') as f:
        json.dump(trimming_params, f, indent=2)
    
    # Print summary
    print(f"Analyzed {len(samples_data)} samples.")
    custom_count = sum(1 for params in trimming_params.values() if params.get('customized', False))
    print(f"Generated customized parameters for {custom_count} samples.")
    print(f"Parameters written to {args.output}")

if __name__ == "__main__":
    main()