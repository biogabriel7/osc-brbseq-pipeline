from snakemake.utils import min_version
import os
import csv

# Ensure minimum Snakemake version for compatibility
min_version("6.0")

# Create essential directories if they don't exist
for dir_path in [
    "logs/fastqc", "logs/multiqc", "logs/trimming", 
    "Analysis/QC/FastQC", "Analysis/QC/MultiQC", "Analysis/QC/Trimming",
    "Analysis/QC/Trimming/Reports", "Analysis/QC/Trimming/MultiQC", 
    "Analysis/Trimmed"
]:
    os.makedirs(dir_path, exist_ok=True)

# Load configuration file (using relative path for portability)
configfile: "resources/config/params.yaml"

# Function to parse metasheet and extract sample information
# This centralizes sample handling for all workflow components
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

# Get samples information
SAMPLES = get_samples_from_metasheet()

# Include rules from other Snakefiles
# This brings in all the rules defined in the other files
include: "workflow/qc_params.snakefile"     # QC and parameter generation rules
include: "workflow/fastp_trimming.snakefile"  # Trimming rules

# Define the workflow order
ruleorder: fastqc > multiqc > generate_trimming_params > fastp_with_dependency > multiqc_fastp

# Define the complete workflow with all expected outputs
rule all:
    input:
        # QC stage outputs
        expand("Analysis/QC/FastQC/{sample}/{sample}_R1_fastqc.html",
               sample=SAMPLES.keys()),
        expand("Analysis/QC/FastQC/{sample}/{sample}_R2_fastqc.html",
               sample=SAMPLES.keys()),
        "Analysis/QC/MultiQC/multiqc_report.html",
        "Analysis/QC/Trimming/trimming_params.json",
        # QC completion marker
        "Analysis/QC/.qc_complete",
        # Trimming stage outputs
        expand("Analysis/Trimmed/{sample}/{sample}_R1_trimmed.fastq.gz",
               sample=SAMPLES.keys()),
        expand("Analysis/Trimmed/{sample}/{sample}_R2_trimmed.fastq.gz",
               sample=SAMPLES.keys()),
        "Analysis/QC/Trimming/MultiQC/multiqc_report.html"

# This rule ensures QC is completed before trimming starts
# It acts as a checkpoint between workflow stages
rule qc_complete:
    input:
        fastqc=expand("Analysis/QC/FastQC/{sample}/{sample}_{read}_fastqc.html",
               sample=SAMPLES.keys(), 
               read=["R1", "R2"]),
        multiqc="Analysis/QC/MultiQC/multiqc_report.html",
        params="Analysis/QC/Trimming/trimming_params.json"
    output:
        touch("Analysis/QC/.qc_complete")
    log:
        "logs/workflow/qc_complete.log"
    shell:
        """
        # Create logs directory if it doesn't exist
        mkdir -p $(dirname {log})
        
        # Log QC completion
        echo "QC and parameter generation completed successfully at $(date)" | tee {log}
        echo "FastQC reports generated: $(find Analysis/QC/FastQC -name '*_fastqc.html' | wc -l)" | tee -a {log}
        echo "MultiQC report: {input.multiqc}" | tee -a {log}
        echo "Trimming parameters: {input.params}" | tee -a {log}
        
        # Verify that all required files exist
        if [ ! -f "{input.multiqc}" ]; then
            echo "ERROR: MultiQC report not found" | tee -a {log}
            exit 1
        fi
        
        if [ ! -f "{input.params}" ]; then
            echo "ERROR: Trimming parameters file not found" | tee -a {log}
            exit 1
        fi
        
        # Create a marker file to indicate QC completion
        echo "Creating QC completion marker file" | tee -a {log}
        """

# Modify the original fastp rule to require QC completion
# This ensures proper stage dependencies
ruleorder: fastp_with_dependency > fastp

rule fastp_with_dependency:
    input:
        r1=lambda wildcards: SAMPLES[wildcards.sample]["R1"],
        r2=lambda wildcards: SAMPLES[wildcards.sample]["R2"],
        qc_complete="Analysis/QC/.qc_complete"
    output:
        r1="Analysis/Trimmed/{sample}/{sample}_R1_trimmed.fastq.gz",
        r2="Analysis/Trimmed/{sample}/{sample}_R2_trimmed.fastq.gz",
        html="Analysis/QC/Trimming/Reports/{sample}_fastp.html",
        json="Analysis/QC/Trimming/Reports/{sample}_fastp.json"
    log:
        "logs/trimming/{sample}.log"
    params:
        # Use only the sample name
        sample_params=lambda wildcards: get_sample_params(wildcards.sample),
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
        
        # Extract parameters
        PARAMS='{params.sample_params}'
        echo "Parameters: $PARAMS" >> {log}
        
        # Extract quality settings using a temporary JSON file to avoid shell quoting issues
        echo $PARAMS > params_temp.json
        LEADING_QUALITY=$(python -c "import json; print(json.load(open('params_temp.json')).get('leading_quality', 3))")
        TRAILING_QUALITY=$(python -c "import json; print(json.load(open('params_temp.json')).get('trailing_quality', 3))")
        MIN_LENGTH=$(python -c "import json; print(json.load(open('params_temp.json')).get('min_length', 36))")
        
        # Parse sliding window into components
        SLIDING_WINDOW=$(python -c "import json; print(json.load(open('params_temp.json')).get('sliding_window', '4:15'))")
        WINDOW_SIZE=$(echo $SLIDING_WINDOW | cut -d: -f1)
        WINDOW_QUALITY=$(echo $SLIDING_WINDOW | cut -d: -f2)
        
        echo "Extracted parameters:" >> {log}
        echo "  Leading quality: $LEADING_QUALITY" >> {log}
        echo "  Trailing quality: $TRAILING_QUALITY" >> {log}
        echo "  Min length: $MIN_LENGTH" >> {log}
        echo "  Window size: $WINDOW_SIZE" >> {log}
        echo "  Window quality: $WINDOW_QUALITY" >> {log}
        
        # Build custom options
        CUSTOM_OPTS=""
        
        # Check for specific flags
        DEDUP=$(python -c "import json; print('true' if json.load(open('params_temp.json')).get('deduplication', False) else 'false')")
        if [ "$DEDUP" = "true" ]; then
            CUSTOM_OPTS="$CUSTOM_OPTS --dedup"
            echo "  Using deduplication" >> {log}
        fi
        
        TRIM_POLY_G=$(python -c "import json; print('true' if json.load(open('params_temp.json')).get('trim_poly_g', False) else 'false')")
        if [ "$TRIM_POLY_G" = "true" ]; then
            CUSTOM_OPTS="$CUSTOM_OPTS --trim_poly_g"
            echo "  Using poly-G trimming" >> {log}
        fi
        
        LOW_COMPLEXITY=$(python -c "import json; print('true' if json.load(open('params_temp.json')).get('low_complexity_filter', False) else 'false')")
        if [ "$LOW_COMPLEXITY" = "true" ]; then
            CUSTOM_OPTS="$CUSTOM_OPTS --low_complexity_filter"
            echo "  Using low complexity filter" >> {log}
        fi
        
        # Clean up temporary file
        rm params_temp.json
        
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

# Rules for running specific workflow stages independently
# These allow users to run just a portion of the workflow

# Run just the parameters generation step
rule params_generation:
    input:
        "Analysis/QC/Trimming/trimming_params.json"
    shell:
        """
        echo "Parameter generation completed successfully"
        """

# Run just the QC stage
rule qc_stage:
    input:
        "Analysis/QC/.qc_complete"
    shell:
        """
        echo "QC stage completed successfully"
        """

# Run just the trimming stage
rule trimming_stage:
    input:
        expand("Analysis/Trimmed/{sample}/{sample}_R1_trimmed.fastq.gz",
               sample=SAMPLES.keys()),
        expand("Analysis/Trimmed/{sample}/{sample}_R2_trimmed.fastq.gz",
               sample=SAMPLES.keys()),
        "Analysis/QC/Trimming/MultiQC/multiqc_report.html"
    shell:
        """
        echo "Trimming stage completed successfully"
        """

# Define workflow stages for better control over execution order
rule workflow_stages:
    input:
        qc="Analysis/QC/.qc_complete",
        trimming=expand("Analysis/Trimmed/{sample}/{sample}_R1_trimmed.fastq.gz",
                      sample=SAMPLES.keys())
    output:
        touch("Analysis/.workflow_complete")
    shell:
        """
        echo "Workflow completed successfully at $(date)"
        echo "QC stage: Completed"
        echo "Trimming stage: Completed"
        """

# Provide helpful messages on workflow completion or failure
onsuccess:
    print("\nWorkflow completed successfully!")
    print("================================")
    print("Quality control reports: Analysis/QC/MultiQC/multiqc_report.html")
    print("Trimming parameters: Analysis/QC/Trimming/trimming_params.json")
    print("Trimmed reads: Analysis/Trimmed/{sample}/{sample}_R1_trimmed.fastq.gz")
    print("Trimming reports: Analysis/QC/Trimming/Reports/{sample}_fastp.html")
    print("Trimming MultiQC: Analysis/QC/Trimming/MultiQC/multiqc_report.html")
    print("\nTo visualize the workflow graph:")
    print("  snakemake --dag | dot -Tsvg > workflow.svg")

onerror:
    print("\nWorkflow failed!")
    print("===============")
    print("Check log files in the logs/ directory for error details.")
    print("Run with --notemp flag to keep temporary files for debugging.")
    print("Use -n (dry run) and -p (print shell commands) for troubleshooting.")