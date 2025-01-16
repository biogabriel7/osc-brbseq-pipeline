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

# Define the final output rule (all samples trimmed)
rule all:
    input:
        expand("Trimmed/{sample}/{sample}_R1_paired.fastq.gz", sample=SAMPLES.keys()),
        expand("Trimmed/{sample}/{sample}_R2_paired.fastq.gz", sample=SAMPLES.keys())

# Rule to trim the reads for each sample
rule trim_reads:
    input:
        R1=lambda wildcards: SAMPLES[wildcards.sample]["R1"],
        R2=lambda wildcards: SAMPLES[wildcards.sample]["R2"]
    output:
        R1_paired="Trimmed/{sample}/{sample}_R1_paired.fastq.gz",
        R1_unpaired="Trimmed/{sample}/{sample}_R1_unpaired.fastq.gz",
        R2_paired="Trimmed/{sample}/{sample}_R2_paired.fastq.gz",
        R2_unpaired="Trimmed/{sample}/{sample}_R2_unpaired.fastq.gz"
    params:
        threads=config["threads"],
        adapters=config["adapters"],
        clip_settings=config["clip_settings"],
        leading_quality=config["leading_quality"],
        trailing_quality=config["trailing_quality"],
        sliding_window=config["sliding_window"],
        min_length=config["min_length"]
    shell:
        """
        mkdir -p Trimmed/{wildcards.sample} && \
        trimmomatic PE -threads {params.threads} \
            {input.R1} {input.R2} \
            {output.R1_paired} {output.R1_unpaired} \
            {output.R2_paired} {output.R2_unpaired} \
            ILLUMINACLIP:{params.adapters}:{params.clip_settings} \
            LEADING:{params.leading_quality} TRAILING:{params.trailing_quality} \
            SLIDINGWINDOW:{params.sliding_window} MINLEN:{params.min_length} \
            2> Trimmed/{wildcards.sample}/{wildcards.sample}_trimmomatic.log
        """
