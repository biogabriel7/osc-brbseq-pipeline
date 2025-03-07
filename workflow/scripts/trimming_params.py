#!/usr/bin/env python3
"""
Script to analyze MultiQC FastQC results and generate sample-specific trimming parameters.
This script reads multiqc_fastqc.txt, analyzes quality metrics for each sample,
and generates a JSON file with sample-specific trimming parameters.
Optimized for fastp rather than Trimmomatic.

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
import sys
import re
from collections import defaultdict

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Generate trimming parameters from MultiQC FastQC results')
    parser.add_argument('--input', required=True, help='Input multiqc_fastqc.txt file')
    parser.add_argument('--output', required=True, help='Output JSON file for trimming parameters')
    parser.add_argument('--config', help='Snakemake config YAML file with default parameters')
    parser.add_argument('--debug', action='store_true', help='Enable verbose debug output')
    return parser.parse_args()

def setup_logging(debug=False):
    """Set up basic logging."""
    import logging
    level = logging.DEBUG if debug else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    return logging.getLogger('trimming_params')

def load_default_params(config_file, logger):
    """Load default trimming parameters from config file."""
    if not config_file or not os.path.exists(config_file):
        logger.warning(f"No config file found at {config_file}. Using built-in defaults.")
        # Default parameters if config file not provided
        return {
            "leading_quality": 3,
            "trailing_quality": 3,
            "sliding_window": "4:15",
            "min_length": 36
        }
    
    logger.info(f"Loading parameters from config file: {config_file}")
    
    try:
        with open(config_file, 'r') as f:
            config = yaml.safe_load(f)
        
        # Extract trimming parameters from config
        params = {
            "leading_quality": config.get("leading_quality", 3),
            "trailing_quality": config.get("trailing_quality", 3),
            "sliding_window": config.get("sliding_window", "4:15"),
            "min_length": config.get("min_length", 36)
        }
        
        logger.info(f"Loaded default parameters: {params}")
        return params
    
    except Exception as e:
        logger.error(f"Error loading config file: {str(e)}")
        # Fall back to defaults
        logger.warning("Using built-in default parameters")
        return {
            "leading_quality": 3,
            "trailing_quality": 3,
            "sliding_window": "4:15",
            "min_length": 36
        }

def read_fastqc_data(input_file, logger):
    """Read and parse the multiqc_fastqc.txt file."""
    samples_data = {}
    
    logger.info(f"Reading FastQC data from: {input_file}")
    
    try:
        with open(input_file, 'r') as f:
            # First, check file format
            header_line = f.readline().strip()
            f.seek(0)  # Reset to beginning of file
            
            if not header_line or '\t' not in header_line:
                logger.error(f"File appears to be malformed. First line: '{header_line}'")
                raise ValueError("Input file is not in the expected tab-delimited format")
            
            # Try to determine the column names
            header_columns = header_line.split('\t')
            logger.info(f"Detected columns: {header_columns}")
            
            sample_column = None
            if 'Sample' in header_columns:
                sample_column = 'Sample'
            elif 'Filename' in header_columns:
                sample_column = 'Filename'
            else:
                logger.error("Could not find Sample or Filename column in the header")
                raise ValueError("Missing required column in input file")
            
            logger.info(f"Using '{sample_column}' column to identify samples")
            
            # Read the data
            reader = csv.DictReader(f, delimiter='\t')
            
            # Process each row
            for i, row in enumerate(reader):
                if sample_column not in row:
                    logger.warning(f"Row {i+1} missing '{sample_column}' column. Skipping.")
                    continue
                
                # Extract sample name
                sample_name = row[sample_column]
                logger.debug(f"Processing row for: {sample_name}")
                
                # Extract read direction (R1 or R2)
                if '_R1' in sample_name:
                    read_direction = 'R1'
                elif '_R2' in sample_name:
                    read_direction = 'R2'
                else:
                    read_direction = 'unknown'
                    logger.warning(f"Could not determine read direction for: {sample_name}")
                
                # Extract SRR ID or base sample name
                base_sample = None
                
                # Try to find SRR ID pattern
                srr_match = re.search(r'(SRR\d+)', sample_name)
                if srr_match:
                    base_sample = srr_match.group(1)
                    logger.debug(f"Found SRR ID: {base_sample} in {sample_name}")
                else:
                    # Extract base sample name by removing read direction
                    base_sample = sample_name
                    if '_R1' in sample_name:
                        base_sample = sample_name.replace('_R1', '')
                    elif '_R2' in sample_name:
                        base_sample = sample_name.replace('_R2', '')
                    
                    # Try to extract from filename pattern
                    filename_match = re.search(r'(.+?)_R[12]', sample_name)
                    if filename_match:
                        base_sample = filename_match.group(1)
                    
                    logger.debug(f"Extracted base sample name: {base_sample} from {sample_name}")
                
                # Initialize sample data if not present
                if base_sample not in samples_data:
                    samples_data[base_sample] = {}
                
                # Store the row data with original case preserved
                samples_data[base_sample][read_direction] = row
                logger.debug(f"Stored data for {base_sample}/{read_direction}")
            
            # Summary information
            logger.info(f"Processed {len(samples_data)} samples from {input_file}")
            for sample in samples_data:
                read_dirs = list(samples_data[sample].keys())
                logger.info(f"Sample: {sample}, Read directions: {read_dirs}")
            
            if not samples_data:
                logger.warning("No samples were found in the input file")
                # Print file preview for debugging
                with open(input_file, 'r') as f:
                    logger.warning("File content preview (first 5 lines):")
                    for i, line in enumerate(f):
                        if i < 5:
                            logger.warning(f"Line {i+1}: {line.strip()}")
                        else:
                            break
    
    except Exception as e:
        logger.error(f"Error reading FastQC data: {str(e)}")
        raise
    
    return samples_data

def analyze_sample(sample_name, sample_data, default_params, logger):
    """
    Analyze FastQC metrics for a sample and determine appropriate trimming parameters.
    
    Args:
        sample_name: Name of the sample
        sample_data: Dict containing R1 and R2 FastQC metrics
        default_params: Dict containing default trimming parameters
        logger: Logger object
    
    Returns:
        Dict with sample-specific trimming parameters
    """
    # Start with default parameters
    params = default_params.copy()
    
    logger.info(f"Analyzing sample: {sample_name}")
    
    # Get read data
    r1_data = sample_data.get('R1', {})
    r2_data = sample_data.get('R2', {})
    
    # Debug output
    logger.debug(f"R1 data keys: {list(r1_data.keys())}")
    logger.debug(f"R2 data keys: {list(r2_data.keys())}")
    
    # Flag to track if any customization was applied
    customized = False
    
    # Helper function to find case-insensitive key in dictionary
    def find_key(data, key_name):
        return next((k for k in data.keys() if k.lower() == key_name.lower()), None)
    
    # Check for adapter content issues
    adapter_key = find_key(r1_data, 'adapter_content')
    if adapter_key:
        adapter_r1 = r1_data.get(adapter_key)
        adapter_r2 = r2_data.get(adapter_key) if adapter_key in r2_data else None
        
        logger.debug(f"Adapter content - R1: {adapter_r1}, R2: {adapter_r2}")
        
        if adapter_r1 in ['warn', 'fail'] or adapter_r2 in ['warn', 'fail']:
            # With fastp, we're using auto-detection, but we could adjust stringency if needed
            params['adapter_stringency'] = 0.8  # For reference - fastp handles this differently
            logger.info(f"Detected adapter issues, setting stringency to 0.8")
            customized = True
    
    # Check for quality issues
    quality_key = find_key(r1_data, 'per_base_sequence_quality')
    if quality_key:
        quality_r1 = r1_data.get(quality_key)
        quality_r2 = r2_data.get(quality_key) if quality_key in r2_data else None
        
        logger.debug(f"Sequence quality - R1: {quality_r1}, R2: {quality_r2}")
        
        if quality_r1 in ['warn', 'fail'] or quality_r2 in ['warn', 'fail']:
            # More aggressive quality trimming
            params['leading_quality'] = max(5, params['leading_quality'])
            params['trailing_quality'] = max(5, params['trailing_quality'])
            params['sliding_window'] = "4:20"  # Higher quality threshold
            logger.info(f"Detected quality issues, increasing quality thresholds")
            customized = True
    
    # Check for sequence length distribution issues
    length_key = find_key(r1_data, 'sequence_length_distribution')
    if length_key:
        length_r1 = r1_data.get(length_key)
        length_r2 = r2_data.get(length_key) if length_key in r2_data else None
        
        logger.debug(f"Length distribution - R1: {length_r1}, R2: {length_r2}")
        
        if length_r1 in ['warn', 'fail'] or length_r2 in ['warn', 'fail']:
            # Adjust minimum length requirement
            seq_length_key = find_key(r1_data, 'avg_sequence_length')
            
            if seq_length_key and r1_data.get(seq_length_key):
                try:
                    # Try to convert to float first
                    seq_length = float(r1_data.get(seq_length_key, 150))
                    # Then convert to int for calculation
                    params['min_length'] = max(36, int(seq_length * 0.7))  # At least 70% of average length
                    logger.info(f"Setting min_length to {params['min_length']} based on avg length {seq_length}")
                except (ValueError, TypeError):
                    logger.warning(f"Could not convert sequence length '{r1_data.get(seq_length_key)}' to number")
                    params['min_length'] = 36
            else:
                # If no length info, use default
                params['min_length'] = 36
                
            customized = True
    
    # Check for N content
    n_content_key = find_key(r1_data, 'per_base_n_content')
    if n_content_key:
        n_content_r1 = r1_data.get(n_content_key)
        n_content_r2 = r2_data.get(n_content_key) if n_content_key in r2_data else None
        
        logger.debug(f"N content - R1: {n_content_r1}, R2: {n_content_r2}")
        
        if n_content_r1 in ['warn', 'fail'] or n_content_r2 in ['warn', 'fail']:
            # Max N filter - fastp has n_base_limit parameter
            params['n_base_limit'] = 5  # Limit bases with N
            logger.info(f"Detected high N content, setting n_base_limit to 5")
            customized = True
    
    # Check for duplication levels
    duplication_key = find_key(r1_data, 'sequence_duplication_levels')
    dedup_pct_key = find_key(r1_data, 'total_deduplicated_percentage')
    
    if duplication_key and dedup_pct_key:
        duplication_r1 = r1_data.get(duplication_key)
        
        logger.debug(f"Duplication levels - R1: {duplication_r1}")
        
        if duplication_r1 in ['warn', 'fail']:
            # Check deduplication percentage
            try:
                total_dedup_pct = float(r1_data.get(dedup_pct_key, 100))
                logger.debug(f"Deduplication percentage: {total_dedup_pct}%")
                
                if total_dedup_pct < 50:  # High duplication rate
                    params['deduplication'] = True
                    logger.info(f"High duplication detected ({total_dedup_pct}%), enabling deduplication")
                    customized = True
            except (ValueError, TypeError):
                logger.warning(f"Could not convert deduplication percentage '{r1_data.get(dedup_pct_key)}' to number")
    
    # Check for overrepresented sequences
    overrep_key = find_key(r1_data, 'overrepresented_sequences')
    if overrep_key:
        overrep_r1 = r1_data.get(overrep_key)
        overrep_r2 = r2_data.get(overrep_key) if overrep_key in r2_data else None
        
        logger.debug(f"Overrepresented sequences - R1: {overrep_r1}, R2: {overrep_r2}")
        
        if overrep_r1 in ['warn', 'fail'] or overrep_r2 in ['warn', 'fail']:
            # Use fastp's low complexity filter
            params['low_complexity_filter'] = True
            logger.info(f"Detected overrepresented sequences, enabling low complexity filter")
            customized = True
    
    # Check for per base sequence content issues
    base_content_key = find_key(r1_data, 'per_base_sequence_content')
    if base_content_key:
        base_content_r1 = r1_data.get(base_content_key)
        base_content_r2 = r2_data.get(base_content_key) if base_content_key in r2_data else None
        
        logger.debug(f"Base content - R1: {base_content_r1}, R2: {base_content_r2}")
        
        if base_content_r1 in ['warn', 'fail'] or base_content_r2 in ['warn', 'fail']:
            # This could be due to various factors
            if "NextSeq" in sample_name or "NovaSeq" in sample_name:
                params['trim_poly_g'] = True
                logger.info(f"Sample name suggests NextSeq/NovaSeq, enabling poly-G trimming")
                customized = True
                
            # Check if GC content is very skewed
            gc_key = find_key(r1_data, '%gc')
            if gc_key:
                try:
                    gc_content = float(r1_data.get(gc_key, 50))
                    logger.debug(f"GC content: {gc_content}%")
                    
                    if gc_content > 60 or gc_content < 40:
                        params['correction'] = True  # Enable correction
                        logger.info(f"Skewed GC content ({gc_content}%), enabling correction")
                        customized = True
                except (ValueError, TypeError):
                    logger.warning(f"Could not convert GC content '{r1_data.get(gc_key)}' to number")
    
    # Tag whether this sample needed customization
    params['customized'] = customized
    
    if customized:
        logger.info(f"Generated custom parameters for {sample_name}")
    else:
        logger.info(f"Using default parameters for {sample_name}")
    
    return params

def generate_trimming_params(samples_data, default_params, logger):
    """Generate sample-specific trimming parameters based on QC metrics."""
    trimming_params = {}
    
    logger.info(f"Generating parameters for {len(samples_data)} samples")
    
    for sample_name, sample_data in samples_data.items():
        # Analyze sample and get appropriate parameters
        params = analyze_sample(sample_name, sample_data, default_params, logger)
        
        # Store parameters for this sample
        trimming_params[sample_name] = params
    
    return trimming_params

def main():
    """Main function to run the script."""
    args = parse_args()
    
    # Set up logging
    logger = setup_logging(args.debug)
    
    logger.info("Starting trimming parameter generation")
    logger.info(f"Input file: {args.input}")
    logger.info(f"Output file: {args.output}")
    logger.info(f"Config file: {args.config or 'Not provided, using defaults'}")
    
    # Load default parameters
    default_params = load_default_params(args.config, logger)
    
    # Verify input file exists
    if not os.path.exists(args.input):
        logger.error(f"Input file does not exist: {args.input}")
        sys.exit(1)
    
    # Read FastQC data
    try:
        samples_data = read_fastqc_data(args.input, logger)
    except Exception as e:
        logger.error(f"Failed to read FastQC data: {str(e)}")
        sys.exit(1)
    
    # Generate trimming parameters
    trimming_params = generate_trimming_params(samples_data, default_params, logger)
    
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        logger.info(f"Creating output directory: {output_dir}")
        os.makedirs(output_dir, exist_ok=True)
    
    # Output JSON file
    try:
        with open(args.output, 'w') as f:
            json.dump(trimming_params, f, indent=2)
        logger.info(f"Parameters written to {args.output}")
    except Exception as e:
        logger.error(f"Failed to write output file: {str(e)}")
        sys.exit(1)
    
    # Print summary
    custom_count = sum(1 for params in trimming_params.values() if params.get('customized', False))
    logger.info(f"Analyzed {len(samples_data)} samples")
    logger.info(f"Generated customized parameters for {custom_count} samples")
    logger.info("Parameter generation complete")

if __name__ == "__main__":
    main()