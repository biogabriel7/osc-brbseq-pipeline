import os
from snakemake.utils import min_version
min_version("6.0")

# Note: No configfile statement here anymore - config is loaded in main Snakefile

# Rule to generate all QC outputs
rule qc_all:
    input:
        # FastQC outputs
        expand("Analysis/QC/FastQC/{sample}/{sample}_{read}_fastqc.html",
               sample=SAMPLES.keys(), 
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
        outdir="Analysis/QC/FastQC/{sample}"
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
        outdir="Analysis/QC/MultiQC"
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
rule generate_trimming_params:
    input:
        multiqc_html="Analysis/QC/MultiQC/multiqc_report.html",
        config=config.get("params_path", "resources/config/params.yaml")
    output:
        params="Analysis/QC/Trimming/trimming_params.json"
    log:
        "logs/trimming/generate_params.log"
    resources:
        mem_mb=config.get("params_memory", 4000)
    shell:
        """
        # Define the path to the MultiQC FastQC metrics file
        MULTIQC_FASTQC_FILE="Analysis/QC/MultiQC/multiqc_data/multiqc_fastqc.txt"
        
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