#!/bin/bash
#SBATCH --job-name=snakemake_rnaseq
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --mem=64GB
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
mkdir -p logs/star
mkdir -p resources/metadata

# Print workflow information
echo "Starting RNA-seq workflow at $(date)"
echo "Working directory: $(pwd)"
echo "Snakefile: $(realpath Snakefile)"
echo "Config file: $(realpath resources/config/params.yaml)"

# Print environment information
echo "Python version:"
python --version
echo "Conda environment:"
conda info
echo "Available disk space:"
df -h .

# Run Snakemake with cluster configuration
snakemake --snakefile Snakefile \
    --cores 28 \
    --jobs 10 \
    --use-conda \
    --latency-wait 600 \
    --restart-times 3 \
    --keep-going \
    --rerun-incomplete \
    --printshellcmds \
    --reason \
    --verbose \
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
    echo "Trimmed samples CSV: resources/metadata/trimmed_samples.csv"
    echo "Alignment BAMs: Analysis/Alignment/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam"
    echo "Alignment counts: Analysis/Alignment/STAR/{sample}/{sample}.ReadsPerGene.out.tab"
    echo "Alignment reports: Analysis/Alignment/MultiQC/multiqc_report.html"
else
    echo "Workflow failed at $(date)"
    echo "Check logs in logs/slurm/ for details"
    
    # Print additional debugging information
    echo "Checking for critical files:"
    echo "Trimmed samples CSV exists: $(test -f resources/metadata/trimmed_samples.csv && echo 'Yes' || echo 'No')"
    if [ -f resources/metadata/trimmed_samples.csv ]; then
        echo "Trimmed samples CSV content:"
        cat resources/metadata/trimmed_samples.csv
    fi
    
    echo "Trimming complete marker exists: $(test -f Analysis/Trimmed/.trimming_complete && echo 'Yes' || echo 'No')"
    echo "QC complete marker exists: $(test -f Analysis/QC/.qc_complete && echo 'Yes' || echo 'No')"
    
    # List all log files for the failed job
    echo "Recent log files:"
    find logs -type f -name "*.log" -mtime -1 | xargs ls -lh
    
    # Check for OOM errors in log files
    echo "Checking for OOM errors in logs:"
    grep -l "Killed" logs/slurm/*.err | xargs -r echo "OOM detected in:"
    
    exit 1
fi