import csv
import yaml
import os

# Load configuration from the YAML file
with open("resources/config/params.yaml") as f:
    config = yaml.safe_load(f)

# Function to parse metasheet.csv and get sample information
def get_samples_from_metasheet():
    samples = {}
    with open("resources/config/metasheet.csv", "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            sample_name = row["sample"]
            samples[sample_name] = {
                "day": sample_name.split('_')[0]  # Extract MicroD8, MicroD15, etc.
            }
    return samples

# Store samples from metasheet in a dictionary
SAMPLES = get_samples_from_metasheet()

# Group samples by day
def get_samples_by_day(day):
    return [sample for sample in SAMPLES.keys() if f"MicroD{day}" in sample]

# Define sample groups
D8_SAMPLES = get_samples_by_day("8")
D15_SAMPLES = get_samples_by_day("15")
D17_SAMPLES = get_samples_by_day("17")

# Rule to collect all final outputs
rule all:
    input:
        "Analysis/counts/MicroD8_raw_counts.txt",
        "Analysis/counts/MicroD15_raw_counts.txt",
        "Analysis/counts/All_samples_raw_counts.txt"

# Create count matrix for D8 samples
rule create_d8_counts:
    input:
        bams = expand("Analysis/STAR/{sample}/Aligned.sortedByCoord.out.bam", 
                     sample=D8_SAMPLES)
    output:
        "Analysis/counts/MicroD8_raw_counts.txt"
    threads: 8
    resources:
        mem_mb = 32000,
        time = "2:00:00"
    log:
        "logs/counts/d8_featurecounts.log"
    shell:
        """
        module load subread
        featureCounts \
            -T {threads} \
            -t exon \
            -g gene_id \
            -a /users/PAS2598/duarte63/star_index/gencode.v44.basic.annotation.gtf \
            -o {output} \
            {input.bams} \
            2> {log}
        """

# Create count matrix for D15 samples
rule create_d15_counts:
    input:
        bams = expand("Analysis/STAR/{sample}/Aligned.sortedByCoord.out.bam", 
                     sample=D15_SAMPLES)
    output:
        "Analysis/counts/MicroD15_raw_counts.txt"
    threads: 8
    resources:
        mem_mb = 32000,
        time = "2:00:00"
    log:
        "logs/counts/d15_featurecounts.log"
    shell:
        """
        module load subread
        featureCounts \
            -T {threads} \
            -t exon \
            -g gene_id \
            -a /users/PAS2598/duarte63/star_index/gencode.v44.basic.annotation.gtf \
            -o {output} \
            {input.bams} \
            2> {log}
        """

# Create count matrix for all samples
rule create_all_counts:
    input:
        bams = expand("Analysis/STAR/{sample}/Aligned.sortedByCoord.out.bam", 
                     sample=SAMPLES.keys())
    output:
        "Analysis/counts/All_samples_raw_counts.txt"
    threads: 8
    resources:
        mem_mb = 32000,
        time = "2:00:00"
    log:
        "logs/counts/all_samples_featurecounts.log"
    shell:
        """
        module load subread
        featureCounts \
            -T {threads} \
            -t exon \
            -g gene_id \
            -a /users/PAS2598/duarte63/star_index/gencode.v44.basic.annotation.gtf \
            -o {output} \
            {input.bams} \
            2> {log}
        """