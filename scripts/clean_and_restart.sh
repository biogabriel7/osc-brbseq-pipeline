#!/bin/bash

# Script to clean up previous analysis results and restart the workflow
# This ensures a clean execution after fixing syntax issues

echo "=== Cleaning up previous analysis results ==="

# Remove analysis directories
echo "Removing Analysis directory..."
rm -rf Analysis/

# Remove logs
echo "Removing logs..."
rm -rf logs/

# Remove any temporary Snakemake files
echo "Removing Snakemake temporary files..."
rm -rf .snakemake/

# Remove any marker files
echo "Removing marker files..."
rm -f *.complete

# Create essential directories
echo "Creating essential directories..."
mkdir -p logs/slurm
mkdir -p logs/workflow
mkdir -p logs/fastqc
mkdir -p logs/multiqc
mkdir -p logs/trimming
mkdir -p logs/star

echo "=== Clean-up complete ==="
echo "You can now submit the job with: sbatch submit_snakemake.sh" 