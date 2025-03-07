#!/bin/bash
#SBATCH --job-name=snakemake_qc_trim
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --mem=32GB
#SBATCH --account=PAS2598
#SBATCH --output=logs/slurm/snakemake_%j.out
#SBATCH --error=logs/slurm/snakemake_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=gabriel.duarte@osumc.edu

# Load modules
module load miniconda3/24.1.2-py310

# Activate your conda environment
source activate local

# Create directories for logs
mkdir -p logs/slurm
mkdir -p logs/workflow
mkdir -p logs/fastqc
mkdir -p logs/multiqc
mkdir -p logs/trimming

# Print workflow information
echo "Starting RNA-seq workflow at $(date)"
echo "Working directory: $(pwd)"
echo "Snakefile: $(realpath Snakefile)"
echo "Config file: $(realpath resources/config/params.yaml)"

# Run Snakemake with cluster configuration
snakemake --snakefile Snakefile \
    --cores 28 \
    --jobs 10 \
    --use-conda \
    --latency-wait 120 \
    --restart-times 3 \
    --keep-going \
    --rerun-incomplete \
    --printshellcmds \
    --reason \
    --cluster "sbatch \
        --parsable \
        --job-name=sm.{rule} \
        --time={params.time} \
        --mem={resources.mem_mb}M \
        --cpus-per-task={threads} \
        --account=PAS2598 \
        --output=logs/slurm/{rule}_%j.out \
        --error=logs/slurm/{rule}_%j.err"

# Check if the workflow completed successfully
if [ $? -eq 0 ]; then
    echo "Workflow completed successfully at $(date)"
    
    # Generate workflow visualization
    snakemake --snakefile Snakefile --dag | dot -Tsvg > workflow_dag.svg
    echo "Workflow visualization generated: workflow_dag.svg"
    
    # Print summary of outputs
    echo "QC reports: Analysis/QC/MultiQC/multiqc_report.html"
    echo "Trimming parameters: Analysis/QC/Trimming/trimming_params.json"
    echo "Trimmed reads: Analysis/Trimmed/{sample}/{sample}_R1_trimmed.fastq.gz"
    echo "Trimming reports: Analysis/QC/Trimming/MultiQC/multiqc_report.html"
else
    echo "Workflow failed at $(date)"
    echo "Check logs in logs/slurm/ for details"
    exit 1
fi