import csv
import yaml
import os

# Load configuration from the YAML file
with open("resources/config/params.yaml") as f:
    config = yaml.safe_load(f)

# Function to parse metasheet.csv and extract sample information
def get_samples_from_metasheet():
    samples = {}
    with open("resources/config/metasheet.csv", "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            sample_name = row["sample"]
            # Now constructing the full paths based on the sample name
            samples[sample_name] = {
                "R1": f"/users/PAS2598/duarte63/bulk_rnaseq/Analysis/Trimmed/{sample_name}/{sample_name}_R1_paired.fastq.gz",
                "R2": f"/users/PAS2598/duarte63/bulk_rnaseq/Analysis/Trimmed/{sample_name}/{sample_name}_R2_paired.fastq.gz"
            }
            print(f"Sample: {sample_name}, R1: {samples[sample_name]['R1']}, R2: {samples[sample_name]['R2']}")
    return samples

# Store samples from metasheet in a dictionary
SAMPLES = get_samples_from_metasheet()

# Rule to collect all final outputs
rule all:
    input:
        # STAR outputs
        expand("Analysis/STAR/{sample}/Aligned.sortedByCoord.out.bam", sample=SAMPLES.keys()),
        expand("Analysis/STAR/{sample}/ReadsPerGene.out.tab", sample=SAMPLES.keys()),
        # Optional: Add alignment stats
        expand("Analysis/STAR/{sample}/alignment_metrics.txt", sample=SAMPLES.keys())

# STAR alignment rule
rule star_align:
    input:
        R1_paired = lambda wildcards: SAMPLES[wildcards.sample]["R1"],
        R2_paired = lambda wildcards: SAMPLES[wildcards.sample]["R2"]
    output:
        bam = "Analysis/STAR/{sample}/Aligned.sortedByCoord.out.bam",
        counts = "Analysis/STAR/{sample}/ReadsPerGene.out.tab",
        stats = "Analysis/STAR/{sample}/alignment_metrics.txt"
    params:
        outdir = "Analysis/STAR/{sample}",
        index = config["star_index"]
    threads: config["threads"]
    resources:
        mem_mb = 120000,
        time = "4:00:00",
        ntasks_per_node = 1
    log:
        "logs/star/{sample}_align.log"
    shell:
        """
        module load star
        mkdir -p {params.outdir}
        
        STAR --runMode alignReads \
             --genomeDir {params.index} \
             --readFilesIn {input.R1_paired} {input.R2_paired} \
             --readFilesCommand zcat \
             --outFileNamePrefix {params.outdir}/ \
             --outSAMtype BAM SortedByCoordinate \
             --outSAMunmapped Within \
             --outSAMattributes Standard \
             --quantMode GeneCounts \
             --runThreadN {threads} \
             --outFilterMultimapNmax 1 \
             --outFilterType BySJout \
             --outSAMstrandField intronMotif \
             2> {log}

        # Index BAM with sambamba
        module load sambamba
        sambamba index -t {threads} {output.bam}
        
        # Generate alignment metrics
        sambamba flagstat -t {threads} {output.bam} > {output.stats}
        """