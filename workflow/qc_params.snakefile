import csv
import yaml
import os
from snakemake.utils import min_version
min_version("6.0")

# Load configuration
configfile: "resources/config/params.yaml"

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
for dir in ["logs/fastqc", "logs/multiqc", "logs/trimming", 
            "Analysis/QC/FastQC", "Analysis/QC/MultiQC", "Analysis/QC/Trimming"]:
    os.makedirs(dir, exist_ok=True)

# Rule to generate all outputs
rule qc_all:
    input:
        # FastQC outputs
        expand("Analysis/QC/FastQC/{sample}/{srr}_{read}_fastqc.html",
               sample=SAMPLES.keys(), 
               srr=[SAMPLES[s]["srr"] for s in SAMPLES.keys()],
               read=["R1", "R2"]),
        # MultiQC report
        "Analysis/QC/MultiQC/multiqc_report.html",
        # Trimming parameters
        "Analysis/QC/Trimming/trimming_params.json"

# Rule to run FastQC on both R1 and R2 for each sample
rule fastqc:
    input:
        r1=lambda wildcards: SAMPLES[wildcards.sample]["R1"],
        r2=lambda wildcards: SAMPLES[wildcards.sample]["R2"]
    output:
        html_r1="Analysis/QC/FastQC/{sample}/{sample}_R1_fastqc.html",
        zip_r1="Analysis/QC/FastQC/{sample}/{sample}_R1_fastqc.zip",
        html_r2="Analysis/QC/FastQC/{sample}/{sample}_R2_fastqc.html",
        zip_r2="Analysis/QC/FastQC/{sample}/{sample}_R2_fastqc.zip"
    log:
        "logs/fastqc/{sample}.log"
    params:
        outdir="Analysis/QC/FastQC/{sample}",
        time=config.get("fastqc_time", "01:00:00")
    threads: config.get("fastqc_threads", 2)
    resources:
        mem_mb=config.get("fastqc_memory", 4000)
    shell:
        """
        set -o pipefail
        mkdir -p {params.outdir}
        
        fastqc --threads {threads} \
            --outdir {params.outdir} \
            {input.r1} {input.r2} \
            > {log} 2>&1
            
        # Rename the files to match expected output filenames
        BASE_R1=$(basename {input.r1} .fastq.gz)
        BASE_R2=$(basename {input.r2} .fastq.gz)
        
        # Move and rename files if needed
        if [ -f "{params.outdir}/${{BASE_R1}}_fastqc.html" ] && [ "{params.outdir}/${{BASE_R1}}_fastqc.html" != "{output.html_r1}" ]; then
            mv "{params.outdir}/${{BASE_R1}}_fastqc.html" "{output.html_r1}"
            mv "{params.outdir}/${{BASE_R1}}_fastqc.zip" "{output.zip_r1}"
        fi
        
        if [ -f "{params.outdir}/${{BASE_R2}}_fastqc.html" ] && [ "{params.outdir}/${{BASE_R2}}_fastqc.html" != "{output.html_r2}" ]; then
            mv "{params.outdir}/${{BASE_R2}}_fastqc.html" "{output.html_r2}"
            mv "{params.outdir}/${{BASE_R2}}_fastqc.zip" "{output.zip_r2}"
        fi
        
        # Ensure the output files exist and have correct permissions
        for f in {output}; do
            if [ ! -f "$f" ]; then
                echo "Expected output file $f was not created" >> {log}
                exit 1
            fi
            chmod 644 "$f"  # Ensure files are readable
        done
        """

# Rule to run MultiQC
rule multiqc:
    input:
        fastqc_outputs=expand("Analysis/QC/FastQC/{sample}/{sample}_{read}_fastqc.zip",
                             sample=SAMPLES.keys(),
                             read=["R1", "R2"])
    output:
        html="Analysis/QC/MultiQC/multiqc_report.html",
        data_dir=directory("Analysis/QC/MultiQC/multiqc_data")
    log:
        "logs/multiqc/multiqc.log"
    params:
        outdir="Analysis/QC/MultiQC",
        time=config.get("multiqc_time", "00:30:00")
    resources:
        mem_mb=config.get("multiqc_memory", 4000)
    shell:
        """
        # Run MultiQC with explicit specification of input files
        multiqc --force \
            --outdir {params.outdir} \
            Analysis/QC/FastQC/ \
            2> {log}
            
        # Ensure the output files exist
        if [ ! -f "{output.html}" ]; then
            echo "ERROR: MultiQC report was not created" >> {log}
            exit 1
        fi
        
        # Also check for the specific metrics file we need
        if [ ! -f "{output.data_dir}/multiqc_fastqc.txt" ]; then
            echo "ERROR: FastQC metrics file was not created by MultiQC" >> {log}
            exit 1
        fi
        """

# Rule to generate trimming parameters
# This rule analyzes FastQC results to determine optimal trimming parameters for each sample
# The resulting JSON file is used by the get_fastp_params function in the main Snakefile
rule generate_trimming_params:
    input:
        multiqc_html="Analysis/QC/MultiQC/multiqc_report.html",
        config=config.get("params_path", "resources/config/params.yaml")
    output:
        params="Analysis/QC/Trimming/trimming_params.json"
    log:
        "logs/trimming/generate_params.log"
    params:
        time=config.get("params_time", "00:30:00")
    shell:
        """
        # Define the path to the MultiQC FastQC metrics file
        MULTIQC_FASTQC_FILE="Analysis/QC/MultiQC/multiqc_data/multiqc_fastqc.txt"
        
        # Add debugging information
        echo "Input multiqc report: {input.multiqc_html}" >> {log}
        echo "MultiQC metrics file: $MULTIQC_FASTQC_FILE" >> {log}
        echo "Input config file: {input.config}" >> {log}
        
        # Make sure the script exists
        if [ ! -f "workflow/scripts/trimming_params.py" ]; then
            echo "ERROR: trimming_params.py script not found" >> {log}
            exit 1
        fi
        
        # Check if the multiqc file exists and has content
        if [ ! -f "$MULTIQC_FASTQC_FILE" ]; then
            echo "ERROR: MultiQC FastQC metrics file not found" >> {log}
            echo "Contents of multiqc_data directory:" >> {log}
            ls -la Analysis/QC/MultiQC/multiqc_data/ >> {log}
            exit 1
        fi
        
        # Display the first few lines of the multiqc file for debugging
        echo "First 10 lines of multiqc file:" >> {log}
        head -n 10 "$MULTIQC_FASTQC_FILE" >> {log}
        
        # Run the parameters generation script with verbose output
        python workflow/scripts/trimming_params.py \
            --input "$MULTIQC_FASTQC_FILE" \
            --output {output.params} \
            --config {input.config} \
            --debug \
            >> {log} 2>&1
            
        # Verify JSON output
        echo "Checking JSON output:" >> {log}
        if [ -f "{output.params}" ]; then
            echo "Output JSON file created successfully" >> {log}
            echo "JSON file size: $(stat -c%s {output.params}) bytes" >> {log}
            # Check if it's valid JSON
            python -c "import json; json.load(open('{output.params}'))" >> {log} 2>&1 || echo "WARNING: Invalid JSON output" >> {log}
        else
            echo "ERROR: JSON output file was not created" >> {log}
            exit 1
        fi
        """