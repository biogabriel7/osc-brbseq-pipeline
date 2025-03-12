#!/usr/bin/env python3
"""
Efficient script to merge featureCounts output files into a single count matrix
and generate normalized counts (CPM).

This script is designed to be called from Snakemake.
"""

import os
import pandas as pd
import numpy as np
import logging
from pathlib import Path

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.FileHandler(snakemake.log[0]), logging.StreamHandler()]
)
logger = logging.getLogger(__name__)

def main():
    """Main function to merge count files."""
    try:
        # Get input and output files from Snakemake
        count_files = snakemake.input.counts
        output_merged = snakemake.output.merged
        output_normalized = snakemake.output.normalized
        
        # Create output directory if it doesn't exist
        os.makedirs(os.path.dirname(output_merged), exist_ok=True)
        
        logger.info(f"Starting to merge {len(count_files)} count files")
        
        # Dictionary to store sample counts
        all_counts = {}
        sample_names = []
        
        # Process each count file
        for count_file in count_files:
            # Extract sample name from directory path
            sample_name = os.path.basename(os.path.dirname(count_file))
            sample_names.append(sample_name)
            
            logger.info(f"Processing sample: {sample_name}")
            
            # Read the count file, skipping the first two header lines
            # featureCounts format: first column is gene ID, 7th column is count
            try:
                df = pd.read_csv(count_file, sep='\t', skiprows=2)
                # Extract gene IDs and counts
                genes = df.iloc[:, 0].values
                counts = df.iloc[:, 6].values
                
                # Store in dictionary
                all_counts[sample_name] = dict(zip(genes, counts))
                
                logger.info(f"Processed {len(genes)} genes for {sample_name}")
            except Exception as e:
                logger.error(f"Error processing {count_file}: {str(e)}")
                raise
        
        # Get all unique gene IDs
        all_genes = sorted(set().union(*[set(counts.keys()) for counts in all_counts.values()]))
        logger.info(f"Total unique genes found: {len(all_genes)}")
        
        # Create merged count matrix
        logger.info("Creating merged count matrix")
        with open(output_merged, 'w') as f:
            # Write header
            f.write("Gene_ID\t" + "\t".join(sample_names) + "\n")
            
            # Write counts for each gene
            for gene in all_genes:
                counts = [str(all_counts[sample].get(gene, 0)) for sample in sample_names]
                f.write(f"{gene}\t" + "\t".join(counts) + "\n")
        
        logger.info(f"Merged count matrix written to {output_merged}")
        
        # Calculate library sizes for normalization
        logger.info("Calculating library sizes for normalization")
        lib_sizes = {}
        for sample in sample_names:
            lib_sizes[sample] = sum(all_counts[sample].values())
            logger.info(f"Library size for {sample}: {lib_sizes[sample]:,}")
        
        # Create normalized count matrix (CPM)
        logger.info("Creating normalized count matrix (CPM)")
        with open(output_normalized, 'w') as f:
            # Write header
            f.write("Gene_ID\t" + "\t".join(sample_names) + "\n")
            
            # Write normalized counts for each gene
            for gene in all_genes:
                # Calculate CPM: (count * 1,000,000) / library_size
                cpms = []
                for sample in sample_names:
                    count = all_counts[sample].get(gene, 0)
                    lib_size = lib_sizes[sample]
                    if lib_size > 0:
                        cpm = (count * 1_000_000) / lib_size
                    else:
                        cpm = 0
                    cpms.append(f"{cpm:.2f}")
                
                f.write(f"{gene}\t" + "\t".join(cpms) + "\n")
        
        logger.info(f"Normalized count matrix written to {output_normalized}")
        logger.info("Count merging completed successfully")
        
    except Exception as e:
        logger.error(f"Error in merge_counts.py: {str(e)}")
        raise

if __name__ == "__main__":
    main() 