#!/bin/bash
#SBATCH --job-name=fastqc_pipeline
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28  # Owens has 28 cores per node
#SBATCH --account=PAS2598     # Your project account
#SBATCH --output=logs/slurm/snakemake_%j.out
#SBATCH --error=logs/slurm/snakemake_%j.err

# Create necessary directories
mkdir -p logs/slurm

# Load required modules
module load miniconda3/24.1.2-py310
source activate local

# Execute snakemake
snakemake -s workflow/fastqc.snakefile \
    --cores 28 \
    --jobs 10 \
    --use-conda \
    --conda-frontend mamba \
    --rerun-incomplete \
    --keep-going \
    --latency-wait 60 \
    --printshellcmds \
    --cluster "sbatch \
        --account=PAS2598 \
        --time={resources.time} \
        --mem={resources.mem_mb} \
        --cpus-per-task={threads} \
        --output=logs/slurm/fastqc_%j.out \
        --error=logs/slurm/fastqc_%j.err"