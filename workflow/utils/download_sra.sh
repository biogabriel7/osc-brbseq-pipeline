#!/bin/bash
#SBATCH --account=PAS2598
#SBATCH --time=12:00:00
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=8
#SBATCH --mem=32G
#SBATCH --job-name=prefetch
#SBATCH --output=prefetch_%j.log

# Base directory
BASE_DIR="/fs/ess/PAS2598/protocol_comparison/motoneurons"

# Load conda
module load miniconda3

# Activate your environment
source activate local

# Create main directory and lab directories
mkdir -p ${BASE_DIR}/{Souza,Scaber,Okano}_Samples

# Initialize metasheet
mkdir -p ${BASE_DIR}/resources/config
echo "sample,R1,R2" > ${BASE_DIR}/resources/config/metasheet.csv

# Define samples by lab
declare -A lab_samples=(
    ["Souza"]="SRR22522185 SRR22522186"
    ["Scaber"]="SRR28516486 SRR28516488 SRR28516489"
    ["Okano"]="SRR29476310 SRR29476289 SRR29476277"
)

for lab in "${!lab_samples[@]}"; do
    echo "Processing ${lab} samples..."
    cd ${BASE_DIR}/${lab}_Samples
    
    # Process each sample for this lab
    for srr in ${lab_samples[$lab]}; do
        echo "Prefetching ${srr}..."
        prefetch ${srr}
        
        if [ $? -eq 0 ]; then
            echo "Running fasterq-dump for ${srr}..."
            fasterq-dump --split-files ${srr}
            
            if [ -f "${srr}_1.fastq" ] && [ -f "${srr}_2.fastq" ]; then
                echo "Compressing files..."
                gzip ${srr}_1.fastq
                gzip ${srr}_2.fastq
                
                # Rename files to match desired format
                mv ${srr}_1.fastq.gz ${srr}_R1.fastq.gz
                mv ${srr}_2.fastq.gz ${srr}_R2.fastq.gz
                
                # Add to metasheet
                echo "${lab}_${srr},${BASE_DIR}/${lab}_Samples/${srr}_R1.fastq.gz,${BASE_DIR}/${lab}_Samples/${srr}_R2.fastq.gz" >> ${BASE_DIR}/resources/config/metasheet.csv
                
                # Clean up SRA files
                rm -rf ~/ncbi/public/sra/${srr}.sra
                
                echo "Successfully processed ${srr}"
            else
                echo "Error: FastQ files not created for ${srr}"
            fi
        else
            echo "Error: Prefetch failed for ${srr}"
        fi
    done
done

# Generate download report
{
    echo "Download Report"
    echo "==============="
    echo "Date: $(date)"
    echo ""
    
    echo "Download Summary by Lab:"
    echo "----------------------"
    for lab in "${!lab_samples[@]}"; do
        echo "${lab} Samples:"
        if [ -d "${BASE_DIR}/${lab}_Samples" ]; then
            ls -lh ${BASE_DIR}/${lab}_Samples/*.fastq.gz 2>/dev/null || echo "No files downloaded"
            echo "Total size:"
            du -sh ${BASE_DIR}/${lab}_Samples 2>/dev/null || echo "Directory empty"
        fi
        echo ""
    done
    
    echo "Generated Metasheet Content:"
    echo "========================="
    cat ${BASE_DIR}/resources/config/metasheet.csv
    
    echo ""
    echo "Total Space Used:"
    du -sh ${BASE_DIR}
} > ${BASE_DIR}/download_report.txt

# Deactivate conda environment
conda deactivate

echo "
Download process completed! 
Check ${BASE_DIR}/download_report.txt for detailed information."