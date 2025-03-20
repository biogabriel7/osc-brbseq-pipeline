#!/bin/bash

# Clean up previous analysis results, restart the workflow, and re-run the pipeline


#remove the old merged gene counts files
rm -f Analysis/Counts/FeatureCounts/merged_gene_counts.txt Analysis/Counts/FeatureCounts/merged_gene_counts.normalized.txt
echo "Old merged gene counts files removed"

#run the clean_and_restart.sh script
bash clean_and_restart.sh

#run the submit_snakemake.sh script
bash submit_snakemake.sh
echo "Snakemake job submitted"

echo "Pipeline clean and re-run complete"

#check the job status
squeue -u jmorales558
