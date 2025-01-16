import csv
import yaml
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
            samples[sample_name] = {
                "R1": row["R1"],
                "R2": row["R2"]
            }
            print(f"Sample: {sample_name}, R1: {samples[sample_name]['R1']}, R2: {samples[sample_name]['R2']}")
    return samples
# Store samples from metasheet in a dictionary
SAMPLES = get_samples_from_metasheet()
# Rule to collect all final outputs
rule all:
    input:
        expand("Analysis/Aligned/{sample}/{sample}_sorted.bam", sample=SAMPLES.keys())
# Rule to align reads with HISAT2
rule align_reads:
    input:
        R1_paired="Analysis/Trimmed/{sample}/{sample}_R1_paired.fastq.gz",
        R2_paired="Analysis/Trimmed/{sample}/{sample}_R2_paired.fastq.gz"
    output:
        bam="Analysis/Aligned/{sample}/{sample}.bam"
    params:
        index=config["hisat2_index"]
    threads: config["threads"]
    resources:
        mem_mb=16000,  # 16GB memory
        time="04:00:00"  # 4 hours runtime
    shell:
        """
        mkdir -p Analysis/Aligned && \
        hisat2 -p {threads} -x {params.index} \
            -1 {input.R1_paired} -2 {input.R2_paired} | \
            samtools view -bS - > {output.bam}
        """
# Rule to sort BAM files using Sambamba
rule sort_bam:
    input:
        bam="Analysis/Aligned/{sample}/{sample}.bam"
    output:
        sorted_bam="Analysis/Aligned/{sample}/{sample}_sorted.bam"
    threads: config["threads"]
    resources:
        mem_mb=32000,  # 32GB memory
        time="02:00:00"  # 2 hours runtime
    params:
        temp_dir="/tmp"  # Adjust this path to use a fast storage
    shell:
        """
        sambamba sort -t {threads} -m {resources.mem_mb}M \
            --tmpdir={params.temp_dir} \
            -o {output.sorted_bam} {input.bam} && \
        rm {input.bam}
        """