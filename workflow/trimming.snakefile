import json
import yaml
import os
import csv
from snakemake.utils import min_version
min_version("6.0")

BASE_DIR = config.get("base_dir", os.getcwd())

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
    except FileNotFoundError:
        print(f"Error: Metasheet not found at {config['metasheet_path']}")
        raise
    except KeyError as e:
        print(f"Error: Missing required column in metasheet: {e}")
        raise
    return samples

# Store samples from metasheet in a dictionary
SAMPLES = get_samples_from_metasheet()

# Load sample-specific trimming parameters
def get_trimming_params():
    try:
        with open("Analysis/QC/Trimming/trimming_params.json", "r") as f:
            return json.load(f)
    except FileNotFoundError:
        # Return empty dict if file doesn't exist yet
        return {}

# Create the directory structure
for dir in ["logs/trimming", "Analysis/QC/Trimming", "Analysis/Trimmed"]:
    os.makedirs(dir, exist_ok=True)

# Rule to generate trimming parameters
rule generate_trimming_params:
    input:
        # Use absolute path to be certain
        multiqc_fastqc=os.path.join(BASE_DIR, "Analysis/QC/MultiQC/multiqc_data/multiqc_fastqc.txt"),
        config=config["params_path"]
    output:
        params="Analysis/QC/Trimming/trimming_params.json"
    log:
        "logs/trimming/generate_params.log"
    shell:
        """
        mkdir -p workflow/scripts
        python {BASE_DIR}/workflow/scripts/trimming_params.py \
            --input {input.multiqc_fastqc} \
            --output {output.params} \
            --config {input.config} \
            > {log} 2>&1
        """

# Rule to run Trimmomatic with sample-specific parameters
rule trimmomatic:
    input:
        r1=lambda wildcards: SAMPLES[wildcards.sample]["R1"],
        r2=lambda wildcards: SAMPLES[wildcards.sample]["R2"],
        params="Analysis/QC/Trimming/trimming_params.json"
    output:
        r1="Analysis/Trimmed/{sample}/{srr}_R1_trimmed.fastq.gz",
        r1_unpaired="Analysis/Trimmed/{sample}/{srr}_R1_unpaired.fastq.gz",
        r2="Analysis/Trimmed/{sample}/{srr}_R2_trimmed.fastq.gz",
        r2_unpaired="Analysis/Trimmed/{sample}/{srr}_R2_unpaired.fastq.gz"
    log:
        "logs/trimming/{sample}_{srr}.log"
    params:
        # Get sample-specific parameters or default to config values
        get_params=lambda wildcards, input: get_sample_params(wildcards.sample, input.params)
    threads: config.get("trimmomatic_threads", 4)
    resources:
        mem_mb=config.get("trimmomatic_memory", 8000),
        time=config.get("trimmomatic_time", "02:00:00")
    shell:
        """
        # Extract parameters for this sample
        LEADING=$(echo {params.get_params} | jq -r '.leading_quality')
        TRAILING=$(echo {params.get_params} | jq -r '.trailing_quality')
        SLIDINGWINDOW=$(echo {params.get_params} | jq -r '.sliding_window')
        MINLEN=$(echo {params.get_params} | jq -r '.min_length')
        ADAPTERS=$(echo {params.get_params} | jq -r '.adapters')
        
        # Check if we need to use default ILLUMINACLIP or custom settings
        if [ "$(echo {params.get_params} | jq -r '.customized // false')" = "true" ]; then
            # Customize adapter parameters if needed
            ADAPTER_STRINGENCY=$(echo {params.get_params} | jq -r '.adapter_stringency // 0.7')
            MIN_ADAPTER_LEN=$(echo {params.get_params} | jq -r '.min_adapter_length // 8')
            ADAPTER_CLIP="ILLUMINACLIP:$ADAPTERS:2:$MIN_ADAPTER_LEN:10:$ADAPTER_STRINGENCY"
        else
            # Use default adapter settings
            ADAPTER_CLIP="ILLUMINACLIP:$ADAPTERS:2:30:10:2:True"
        fi
        
        # Add AVGQUAL filter if specific quality issues detected
        if [ "$(echo {params.get_params} | jq -r '.per_base_sequence_quality // "pass"')" != "pass" ]; then
            AVGQUAL="AVGQUAL:20"
        else
            AVGQUAL=""
        fi
        
        # Add MAXINFO parameter for length distribution issues
        if [ "$(echo {params.get_params} | jq -r '.sequence_length_distribution // "pass"')" != "pass" ]; then
            MAXINFO="MAXINFO:40:0.5"
        else
            MAXINFO=""
        fi
                
        # Run Trimmomatic with the determined parameters
        trimmomatic PE -threads {threads} \
            {input.r1} {input.r2} \
            {output.r1} {output.r1_unpaired} \
            {output.r2} {output.r2_unpaired} \
            $ADAPTER_CLIP \
            LEADING:$LEADING \
            TRAILING:$TRAILING \
            SLIDINGWINDOW:$SLIDINGWINDOW \
            $MAXINFO \
            $AVGQUAL \
            MINLEN:$MINLEN \
            > {log} 2>&1
        """

# Helper function to extract sample-specific parameters
def get_sample_params(sample_name, params_file):
    try:
        with open(params_file, "r") as f:
            all_params = json.load(f)
            
        # Get parameters for this sample, or use defaults
        if sample_name in all_params:
            return json.dumps(all_params[sample_name])
        else:
            # Return default parameters from config
            return json.dumps({
                "adapters": config.get("adapters", "resources/adaptors/TruSeq3-PE.fa/TruSeq3-PE.fa"),
                "leading_quality": config.get("leading_quality", 3),
                "trailing_quality": config.get("trailing_quality", 3),
                "sliding_window": config.get("sliding_window", "4:15"),
                "min_length": config.get("min_length", 36),
                "customized": False
            })
    except Exception as e:
        print(f"Error loading parameters for {sample_name}: {e}")
        # Return defaults as fallback
        return json.dumps({
            "adapters": config.get("adapters", "resources/adaptors/TruSeq3-PE.fa/TruSeq3-PE.fa"),
            "leading_quality": config.get("leading_quality", 3),
            "trailing_quality": config.get("trailing_quality", 3),
            "sliding_window": config.get("sliding_window", "4:15"),
            "min_length": config.get("min_length", 36),
            "customized": False
        })

# Rule to run trimming for all samples
rule trim_all:
    input:
        expand("Analysis/Trimmed/{sample}/{srr}_R1_trimmed.fastq.gz",
               sample=SAMPLES.keys(), 
               srr=[SAMPLES[s]["srr"] for s in SAMPLES.keys()]),
        expand("Analysis/Trimmed/{sample}/{srr}_R2_trimmed.fastq.gz",
               sample=SAMPLES.keys(), 
               srr=[SAMPLES[s]["srr"] for s in SAMPLES.keys()])

# Main rule 
rule all:
    input:
        "Analysis/QC/Trimming/trimming_params.json",
        rules.trim_all.input