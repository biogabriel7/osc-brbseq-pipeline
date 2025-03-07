import json
import yaml
import os
import csv
import re
from snakemake.utils import min_version
min_version("6.0")

# Load configuration
configfile: "resources/config/params.yaml"

# Create the directory structure
for dir in ["logs/trimming", "Analysis/Trimmed", "Analysis/QC/Trimming/Reports", "Analysis/QC/Trimming/MultiQC"]:
    os.makedirs(dir, exist_ok=True)

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
                print(f"Loaded sample: {sample_name}, SRR: {samples[sample_name]['srr']}")
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
for dir in ["logs/trimming", "Analysis/Trimmed", "Analysis/QC/Trimming/Reports"]:
    os.makedirs(dir, exist_ok=True)

# Load the trimming parameters from the JSON file
def load_trimming_params():
    params_file = "Analysis/QC/Trimming/trimming_params.json"
    if not os.path.exists(params_file):
        print(f"WARNING: No trimming parameters file found at {params_file}.")
        print("Will use default parameters for all samples.")
        return {}
    
    try:
        with open(params_file, "r") as f:
            return json.load(f)
    except Exception as e:
        print(f"ERROR: Could not parse trimming parameters: {str(e)}")
        return {}

# Store trimming parameters
TRIMMING_PARAMS = load_trimming_params()

# Rule to run trimming for all samples
rule trim_all:
    input:
        expand("Analysis/Trimmed/{sample}/{sample}_R1_trimmed.fastq.gz",
               sample=SAMPLES.keys()),
        expand("Analysis/Trimmed/{sample}/{sample}_R2_trimmed.fastq.gz",
               sample=SAMPLES.keys()),
        "Analysis/QC/Trimming/MultiQC/multiqc_report.html"

# Rule to run fastp with sample-specific parameters
rule fastp:
    input:
        r1=lambda wildcards: SAMPLES[wildcards.sample]["R1"],
        r2=lambda wildcards: SAMPLES[wildcards.sample]["R2"]
    output:
        r1="Analysis/Trimmed/{sample}/{sample}_R1_trimmed.fastq.gz",
        r2="Analysis/Trimmed/{sample}/{sample}_R2_trimmed.fastq.gz",
        html="Analysis/QC/Trimming/Reports/{sample}_fastp.html",
        json="Analysis/QC/Trimming/Reports/{sample}_fastp.json"
    log:
        "logs/trimming/{sample}.log"
    params:
        # Use the get_fastp_params function from the main Snakefile
        fastp_opts=lambda wildcards: get_fastp_params(wildcards),
        time=config.get("trimming_time", "02:00:00")
    threads: config.get("trimming_threads", 4)
    resources:
        mem_mb=config.get("trimming_memory", 8000)
    shell:
        """
        # Create output directories
        mkdir -p $(dirname {output.r1})
        mkdir -p $(dirname {output.html})
        
        # Debug information
        echo "Processing sample: {wildcards.sample}" >> {log}
        echo "Input R1: {input.r1}" >> {log}
        echo "Input R2: {input.r2}" >> {log}
        echo "Output R1: {output.r1}" >> {log}
        echo "Output R2: {output.r2}" >> {log}
        echo "Parameters: {params.fastp_opts}" >> {log}
        
        # Run fastp with auto-adapter detection and determined parameters
        echo "Running fastp..." >> {log}
        fastp \
            --in1 {input.r1} \
            --in2 {input.r2} \
            --out1 {output.r1} \
            --out2 {output.r2} \
            --html {output.html} \
            --json {output.json} \
            --thread {threads} \
            --detect_adapter_for_pe \
            {params.fastp_opts} \
            >> {log} 2>&1
        
        # Verify outputs
        if [ ! -f "{output.r1}" ] || [ ! -f "{output.r2}" ]; then
            echo "ERROR: Output files were not created properly." >> {log}
            echo "Command failed, check log for details." >> {log}
            exit 1
        fi
        
        echo "Trimming completed successfully" >> {log}
        
        # Print basic stats
        echo "Output file sizes:" >> {log}
        ls -lh {output.r1} {output.r2} >> {log}
        """

# Rule to run MultiQC on fastp reports
rule multiqc_fastp:
    input:
        fastp_reports=expand("Analysis/QC/Trimming/Reports/{sample}_fastp.json",
                           sample=SAMPLES.keys())
    output:
        html="Analysis/QC/Trimming/MultiQC/multiqc_report.html",
        data_dir=directory("Analysis/QC/Trimming/MultiQC/multiqc_data")
    log:
        "logs/trimming/multiqc_fastp.log"
    shell:
        """
        # Create output directory
        mkdir -p $(dirname {output.html})
        
        # Run MultiQC on fastp reports
        echo "Running MultiQC on fastp reports..." >> {log}
        
        # List input files for debugging
        echo "Input files:" >> {log}
        for f in {input.fastp_reports}; do
            echo "  $f" >> {log}
            # Check file size
            ls -lh $f >> {log}
        done
        
        multiqc \
            --force \
            --outdir $(dirname {output.html}) \
            --filename $(basename {output.html}) \
            Analysis/QC/Trimming/Reports/ \
            > {log} 2>&1
            
        # Verify output
        if [ ! -f "{output.html}" ]; then
            echo "ERROR: MultiQC report was not created" >> {log}
            exit 1
        fi
        
        echo "MultiQC completed successfully" >> {log}
        """