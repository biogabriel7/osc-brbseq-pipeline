#!/bin/bash
#SBATCH --job-name=fastp_large
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=64G
#SBATCH --account=PAS2598
#SBATCH --output=logs/slurm/fastp_large_%j.out
#SBATCH --error=logs/slurm/fastp_large_%j.err

# Create logs directory
mkdir -p logs/slurm

# Load required modules
module load miniconda3/24.1.2-py310
source activate local

# Execute snakemake for a specific large sample
# Usage: sbatch submit_large_sample.sh SAMPLE_NAME
# Example: sbatch submit_large_sample.sh Scaber_SRR28516488

if [ -z "$1" ]; then
    echo "Error: No sample name provided"
    echo "Usage: sbatch submit_large_sample.sh SAMPLE_NAME"
    exit 1
fi

SAMPLE=$1

echo "Processing large sample: $SAMPLE with 64GB memory"

# Execute snakemake for the specific sample
snakemake -s Snakefile \
    "Analysis/Trimmed/${SAMPLE}/${SAMPLE}_R1_trimmed.fastq.gz" \
    "Analysis/Trimmed/${SAMPLE}/${SAMPLE}_R2_trimmed.fastq.gz" \
    --cores 4 \
    --printshellcmds \
    --latency-wait 120 \
    --rerun-incomplete \
    --resources mem_mb=64000

echo "Job completed at $(date)"
echo "Check logs/trimming/${SAMPLE}.log for details" 