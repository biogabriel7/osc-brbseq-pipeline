#!/bin/bash
#SBATCH --job-name=star_index
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=64GB
#SBATCH --account=PAS2598
#SBATCH --output=logs/slurm/star_index_%j.out
#SBATCH --error=logs/slurm/star_index_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=gabriel.duarte@osumc.edu

# Load modules
module load miniconda3/24.1.2-py310

# Activate your conda environment
source activate local

# Define paths
GENOME_FASTA="resources/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GTF_FILE="resources/genome/Homo_sapiens.GRCh38.90.gtf"
STAR_INDEX="resources/genome/star_index"
THREADS=16
SJDB_OVERHANG=100  # Read length - 1

# Create output directory
mkdir -p $STAR_INDEX
mkdir -p logs/slurm

# Print information
echo "Creating STAR index at $(date)"
echo "Genome FASTA: $GENOME_FASTA"
echo "GTF file: $GTF_FILE"
echo "Output directory: $STAR_INDEX"
echo "Using $THREADS threads"

# Run STAR index generation
STAR --runMode genomeGenerate \
    --runThreadN $THREADS \
    --genomeDir $STAR_INDEX \
    --genomeFastaFiles $GENOME_FASTA \
    --sjdbGTFfile $GTF_FILE \
    --sjdbOverhang $SJDB_OVERHANG

# Check if the index was created successfully
if [ $? -eq 0 ]; then
    echo "STAR index created successfully at $(date)"
    echo "Index files:"
    ls -lh $STAR_INDEX
else
    echo "STAR index creation failed at $(date)"
    exit 1
fi 