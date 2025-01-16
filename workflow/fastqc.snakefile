import csv
import yaml
import os
from snakemake.utils import min_version
min_version("6.0")

# Load configuration
configfile: "/users/PAS2598/duarte63/GitHub/satienzar-brnaseq/resources/config/params.yaml"

# Function to parse metasheet and extract sample information
def get_samples_from_metasheet():
    samples = {}
    try:
        with open(config["metasheet_path"], "r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                sample_name = row["sample"]
                samples[sample_name] = {
                    "R1": row["R1"],
                    "R2": row["R2"],
                    "srr": os.path.basename(row["R1"]).replace("_R1.fastq.gz", "")  # Get SRR ID from filename
                }
                print(f"Sample: {sample_name}, SRR: {samples[sample_name]['srr']}")
    except FileNotFoundError:
        print(f"Error: Metasheet not found at {config['metasheet_path']}")
        raise
    except KeyError as e:
        print(f"Error: Missing required column in metasheet: {e}")
        raise
    return samples

# Store samples from metasheet in a dictionary
SAMPLES = get_samples_from_metasheet()

# Create the directory structure
for dir in ["logs/fastqc", "Analysis/QC/FastQC"]:
    os.makedirs(dir, exist_ok=True)

# Rule to generate all FastQC reports
rule all:
    input:
        expand("Analysis/QC/FastQC/{sample}/{srr}_{read}_fastqc.html",
               sample=SAMPLES.keys(), 
               srr=[SAMPLES[s]["srr"] for s in SAMPLES.keys()],
               read=["R1", "R2"])

# Rule to run FastQC on both R1 and R2 for each sample
rule fastqc:
    input:
        r1=lambda wildcards: SAMPLES[wildcards.sample]["R1"],
        r2=lambda wildcards: SAMPLES[wildcards.sample]["R2"]
    output:
        html_r1="Analysis/QC/FastQC/{sample}/{srr}_R1_fastqc.html",
        zip_r1="Analysis/QC/FastQC/{sample}/{srr}_R1_fastqc.zip",
        html_r2="Analysis/QC/FastQC/{sample}/{srr}_R2_fastqc.html",
        zip_r2="Analysis/QC/FastQC/{sample}/{srr}_R2_fastqc.zip"
    log:
        "logs/fastqc/{sample}_{srr}.log"
    params:
        outdir="Analysis/QC/FastQC/{sample}"
    threads: config["fastqc_threads"] if "fastqc_threads" in config else 2
    resources:
        mem_mb=config["fastqc_memory"] if "fastqc_memory" in config else 4000,
        time=config["fastqc_time"] if "fastqc_time" in config else "01:00:00"
    shell:
        """
        set -o pipefail
        mkdir -p {params.outdir}
        
        fastqc --threads {threads} \
            --outdir {params.outdir} \
            {input.r1} {input.r2} \
            > {log} 2>&1
            
        # Ensure the output files exist
        for f in {output}; do
            if [ ! -f "$f" ]; then
                echo "Expected output file $f was not created" >> {log}
                exit 1
            fi
        done
        """