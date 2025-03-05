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
        html_r1="Analysis/QC/FastQC/{sample}/{srr}_R1_fastqc.html",
        zip_r1="Analysis/QC/FastQC/{sample}/{srr}_R1_fastqc.zip",
        html_r2="Analysis/QC/FastQC/{sample}/{srr}_R2_fastqc.html",
        zip_r2="Analysis/QC/FastQC/{sample}/{srr}_R2_fastqc.zip"
    log:
        "logs/fastqc/{sample}_{srr}.log"
    params:
        outdir="Analysis/QC/FastQC/{sample}"
    threads: config.get("fastqc_threads", 2)
    resources:
        mem_mb=config.get("fastqc_memory", 4000),
        time=config.get("fastqc_time", "01:00:00")
    shell:
        """
        set -o pipefail
        mkdir -p {params.outdir}
        
        fastqc --threads {threads} \
            --outdir {params.outdir} \
            {input.r1} {input.r2} \
            > {log} 2>&1
            
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
        # Use the directory containing all FastQC outputs
        fastqc_outputs=expand("Analysis/QC/FastQC/{sample}/{srr}_{read}_fastqc.zip",
                             sample=SAMPLES.keys(),
                             srr=[SAMPLES[s]["srr"] for s in SAMPLES.keys()],
                             read=["R1", "R2"])
    output:
        html="Analysis/QC/MultiQC/multiqc_report.html",
        data_dir=directory("Analysis/QC/MultiQC/multiqc_data")
    log:
        "logs/multiqc/multiqc.log"
    params:
        outdir="Analysis/QC/MultiQC"
    resources:
        mem_mb=config.get("multiqc_memory", 4000),
        time=config.get("multiqc_time", "00:30:00")
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
rule generate_trimming_params:
    input:
        multiqc_fastqc="Analysis/QC/MultiQC/multiqc_data/multiqc_fastqc.txt",
        config=config.get("params_path", "resources/config/params.yaml")
    output:
        params="Analysis/QC/Trimming/trimming_params.json"
    log:
        "logs/trimming/generate_params.log"
    shell:
        """
        # Add debugging information
        echo "Input multiqc file: {input.multiqc_fastqc}" >> {log}
        echo "Input config file: {input.config}" >> {log}
        
        # Make sure the script exists
        if [ ! -f "workflow/scripts/trimming_params.py" ]; then
            echo "ERROR: trimming_params.py script not found" >> {log}
            exit 1
        fi
        
        # Check if the multiqc file exists and has content
        if [ ! -f "{input.multiqc_fastqc}" ]; then
            echo "ERROR: MultiQC FastQC metrics file not found" >> {log}
            exit 1
        fi
        
        # Display the first few lines of the multiqc file for debugging
        echo "First 10 lines of multiqc file:" >> {log}
        head -n 10 {input.multiqc_fastqc} >> {log}
        
        # Run the parameters generation script with verbose output
        python workflow/scripts/trimming_params.py \
            --input {input.multiqc_fastqc} \
            --output {output.params} \
            --config {input.config} \
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