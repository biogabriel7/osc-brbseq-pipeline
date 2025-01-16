#!/bin/bash
#SBATCH --job-name=multiqc
#SBATCH --time=0:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --account=PAS2598
#SBATCH --output=logs/slurm/multiqc_%j.out
#SBATCH --error=logs/slurm/multiqc_%j.err

# Create logs directory
mkdir -p logs/slurm

# Load required modules
module load miniconda3/24.1.2-py310
source activate local

# Execute snakemake
snakemake -s workflow/multiqc.snakefile \
    --cores 1 \
    --printshellcmds \
    --latency-wait 60