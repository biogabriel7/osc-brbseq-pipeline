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

# Create the directory structure
for dir in ["logs/trimming", "Analysis/QC/Trimming", "Analysis/Trimmed", "Analysis/QC/Trimming/Reports"]:
    os.makedirs(dir, exist_ok=True)

# Rule to generate trimming parameters
rule generate_trimming_params:
    input:
        multiqc_fastqc=os.path.join(BASE_DIR, "Analysis/QC/MultiQC/multiqc_data/multiqc_fastqc.txt"),
        config=config["params_path"]
    output:
        params="Analysis/QC/Trimming/trimming_params.json"
    log:
        "logs/trimming/generate_params.log"
    shell:
        """
        python {BASE_DIR}/workflow/scripts/trimming_params.py \
            --input {input.multiqc_fastqc} \
            --output {output.params} \
            --config {input.config} \
            > {log} 2>&1
        """

# Helper function to extract sample-specific parameters
def get_sample_params(srr_id, params_file):
    try:
        with open(params_file, "r") as f:
            all_params = json.load(f)
            
        # Get parameters for this SRR ID, or use defaults
        if srr_id in all_params:
            return json.dumps(all_params[srr_id])
        else:
            # Try to look for other variations of the sample name
            # (in case the key format changes in the future)
            for key in all_params:
                # Check if the SRR ID is a substring of the key
                if srr_id in key:
                    return json.dumps(all_params[key])
            
            # Return default parameters from config if no match found
            print(f"Warning: No matching parameters found for {srr_id}. Using defaults.")
            return json.dumps({
                "leading_quality": config.get("leading_quality", 3),
                "trailing_quality": config.get("trailing_quality", 3),
                "sliding_window": config.get("sliding_window", "4:15"),
                "min_length": config.get("min_length", 36),
                "customized": False
            })
    except Exception as e:
        print(f"Error loading parameters for {srr_id}: {e}")
        # Return defaults as fallback
        return json.dumps({
            "leading_quality": config.get("leading_quality", 3),
            "trailing_quality": config.get("trailing_quality", 3),
            "sliding_window": config.get("sliding_window", "4:15"),
            "min_length": config.get("min_length", 36),
            "customized": False
        })

# Rule to run fastp with sample-specific parameters and auto-adapter detection
rule fastp:
    input:
        r1=lambda wildcards: SAMPLES[wildcards.sample]["R1"],
        r2=lambda wildcards: SAMPLES[wildcards.sample]["R2"],
        params="Analysis/QC/Trimming/trimming_params.json"
    output:
        r1="Analysis/Trimmed/{sample}/{srr}_R1_trimmed.fastq.gz",
        r2="Analysis/Trimmed/{sample}/{srr}_R2_trimmed.fastq.gz",
        html="Analysis/QC/Trimming/Reports/{sample}_{srr}_fastp.html",
        json="Analysis/QC/Trimming/Reports/{sample}_{srr}_fastp.json"
    log:
        "logs/trimming/{sample}_{srr}.log"
    params:
        # Get sample-specific parameters using SRR ID instead of sample name
        get_params=lambda wildcards, input: get_sample_params(wildcards.srr, input.params)
    threads: config.get("trimming_threads", 4)
    resources:
        mem_mb=config.get("trimming_memory", 8000),
        time=config.get("trimming_time", "02:00:00")
    shell:
        """
        # Extract parameters for this sample
        PARAMS='{params.get_params}'
        
        # Set quality parameters
        LEADING_QUALITY=$(echo $PARAMS | jq -r '.leading_quality')
        TRAILING_QUALITY=$(echo $PARAMS | jq -r '.trailing_quality')
        MIN_LENGTH=$(echo $PARAMS | jq -r '.min_length')
        
        # Parse sliding window - convert Trimmomatic format to fastp
        SLIDING_WINDOW=$(echo $PARAMS | jq -r '.sliding_window')
        WINDOW_SIZE=$(echo $SLIDING_WINDOW | cut -d: -f1)
        WINDOW_QUALITY=$(echo $SLIDING_WINDOW | cut -d: -f2)
        
        # Detect if customized parameters are needed
        CUSTOM_OPTS=""
        
        # Add deduplication if needed
        if [ "$(echo $PARAMS | jq -r '.deduplication // false')" = "true" ]; then
            CUSTOM_OPTS="$CUSTOM_OPTS --dedup"
        fi
        
        # Add poly-G tail trimming for NextSeq/NovaSeq data if appropriate
        if [ "$(echo $PARAMS | jq -r '.trim_poly_g // false')" = "true" ]; then
            CUSTOM_OPTS="$CUSTOM_OPTS --trim_poly_g"
        fi
        
        # Add correction for low complexity if needed
        if [ "$(echo $PARAMS | jq -r '.low_complexity_filter // false')" = "true" ]; then
            CUSTOM_OPTS="$CUSTOM_OPTS --low_complexity_filter"
        fi
        
        # Debug information
        echo "Processing sample: {wildcards.sample}, SRR: {wildcards.srr}" >> {log}
        echo "Parameters: $PARAMS" >> {log}
        echo "Window size: $WINDOW_SIZE, Window quality: $WINDOW_QUALITY" >> {log}
        echo "Custom options: $CUSTOM_OPTS" >> {log}
        
        # Run fastp with auto-adapter detection and determined parameters
        fastp \
            --in1 {input.r1} \
            --in2 {input.r2} \
            --out1 {output.r1} \
            --out2 {output.r2} \
            --html {output.html} \
            --json {output.json} \
            --thread {threads} \
            --qualified_quality_phred $TRAILING_QUALITY \
            --unqualified_percent_limit 40 \
            --cut_front \
            --cut_front_window_size $WINDOW_SIZE \
            --cut_front_mean_quality $WINDOW_QUALITY \
            --cut_tail \
            --cut_tail_window_size $WINDOW_SIZE \
            --cut_tail_mean_quality $WINDOW_QUALITY \
            --length_required $MIN_LENGTH \
            --detect_adapter_for_pe \
            $CUSTOM_OPTS \
            >> {log} 2>&1
            
        # Check if output files were created successfully
        if [ ! -f "{output.r1}" ] || [ ! -f "{output.r2}" ]; then
            echo "ERROR: Output files were not created properly. Check log for details." >> {log}
            exit 1
        fi
        
        echo "Trimming completed successfully for {wildcards.sample}" >> {log}
        """

# Rule to run MultiQC on fastp reports
rule multiqc_fastp:
    input:
        fastp_reports=expand("Analysis/QC/Trimming/Reports/{sample}_{srr}_fastp.json",
                           sample=SAMPLES.keys(), 
                           srr=[SAMPLES[s]["srr"] for s in SAMPLES.keys()])
    output:
        html="Analysis/QC/Trimming/MultiQC/multiqc_report.html",
        data_dir=directory("Analysis/QC/Trimming/MultiQC/multiqc_data")
    log:
        "logs/trimming/multiqc.log"
    shell:
        """
        # Create output directory if it doesn't exist
        mkdir -p $(dirname {output.html})
        
        # Run MultiQC on fastp reports
        multiqc \
            --force \
            --outdir $(dirname {output.html}) \
            --filename $(basename {output.html}) \
            Analysis/QC/Trimming/Reports/ \
            > {log} 2>&1
        """

# Rule to run trimming for all samples
rule trim_all:
    input:
        expand("Analysis/Trimmed/{sample}/{srr}_R1_trimmed.fastq.gz",
               sample=SAMPLES.keys(), 
               srr=[SAMPLES[s]["srr"] for s in SAMPLES.keys()]),
        expand("Analysis/Trimmed/{sample}/{srr}_R2_trimmed.fastq.gz",
               sample=SAMPLES.keys(), 
               srr=[SAMPLES[s]["srr"] for s in SAMPLES.keys()]),
        "Analysis/QC/Trimming/MultiQC/multiqc_report.html"

# Main rule 
rule all:
    input:
        "Analysis/QC/Trimming/trimming_params.json",
        rules.trim_all.input