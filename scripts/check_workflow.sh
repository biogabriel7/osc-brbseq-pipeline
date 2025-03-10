#!/bin/bash

# Load modules
module load miniconda3/24.1.2-py310

# Activate your conda environment
source activate local

# Print workflow information
echo "Checking RNA-seq workflow at $(date)"
echo "Working directory: $(pwd)"
echo "Snakefile: $(realpath Snakefile)"
echo "Config file: $(realpath resources/config/params.yaml)"

# Run Snakemake in dry-run mode
echo "Running dry-run to check workflow..."
snakemake --snakefile Snakefile \
    --cores 1 \
    --dry-run \
    --printshellcmds \
    --reason

# Check if the dry run was successful
if [ $? -eq 0 ]; then
    echo "Dry run completed successfully. Workflow looks good!"
    echo "You can now submit the job with: sbatch submit_snakemake.sh"
else
    echo "Dry run failed. Please check the errors above."
    exit 1
fi 