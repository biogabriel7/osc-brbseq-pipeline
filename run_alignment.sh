#!/bin/bash
#SBATCH --job-name=alignment_only
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --mem=64GB
#SBATCH --account=PAS2598
#SBATCH --output=logs/slurm/alignment_%j.out
#SBATCH --error=logs/slurm/alignment_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=gabriel.duarte@osumc.edu

# Load modules
module load miniconda3/24.1.2-py310

# Activate your conda environment
source activate local

# Create directories for logs
mkdir -p logs/slurm
mkdir -p logs/star

# Print workflow information
echo "Starting alignment stage at $(date)"
echo "Working directory: $(pwd)"
echo "Snakefile: $(realpath Snakefile)"
echo "Config file: $(realpath resources/config/params.yaml)"

# Check if STAR index exists
STAR_INDEX=$(grep "star_index:" resources/config/params.yaml | awk '{print $2}' | tr -d '"')
if [ ! -d "$STAR_INDEX" ]; then
    echo "ERROR: STAR index not found at $STAR_INDEX"
    echo "Please run the create_star_index.sh script first"
    exit 1
fi

# Run Snakemake with cluster configuration, targeting only the alignment stage
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
    --until alignment_stage \
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
    echo "Alignment stage completed successfully at $(date)"
    
    # Print summary of outputs
    echo "Alignment BAMs: Analysis/Alignment/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam"
    echo "Alignment counts: Analysis/Alignment/STAR/{sample}/{sample}.ReadsPerGene.out.tab"
    echo "Alignment reports: Analysis/Alignment/MultiQC/multiqc_report.html"
else
    echo "Alignment stage failed at $(date)"
    echo "Check logs in logs/slurm/ for details"
    exit 1
fi 