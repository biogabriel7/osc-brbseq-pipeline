from snakemake.utils import min_version
import os
import csv

# Ensure minimum Snakemake version for compatibility
min_version("6.0")

# Create essential directories if they don't exist
for dir_path in [
    "logs/fastqc", "logs/multiqc", "logs/trimming", "logs/star",
    "Analysis/QC/FastQC", "Analysis/QC/MultiQC", "Analysis/QC/Trimming",
    "Analysis/QC/Trimming/Reports", "Analysis/QC/Trimming/MultiQC", 
    "Analysis/Trimmed", "Analysis/Alignment/STAR", "Analysis/Alignment/MultiQC"
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
include: "workflow/star_alignment.snakefile"  # Alignment rules

# Define the workflow order
ruleorder: fastqc > multiqc > generate_trimming_params > fastp_with_dependency > fastp > multiqc_fastp > star_align_with_dependency > star_align > index_bam > alignment_metrics > multiqc_alignment

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
        "Analysis/QC/Trimming/MultiQC/multiqc_report.html",
        "Analysis/QC/Trimming/MultiQC/multiqc_report_data",
        # Trimming completion marker and metadata
        "Analysis/Trimmed/.trimming_complete",
        "resources/metadata/trimmed_samples.csv",
        # Alignment stage outputs
        expand("Analysis/Alignment/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
               sample=SAMPLES.keys()),
        expand("Analysis/Alignment/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai",
               sample=SAMPLES.keys()),
        "Analysis/Alignment/MultiQC/multiqc_report.html",
        # Alignment completion marker
        "Analysis/Alignment/.main_alignment_complete",
        # Final workflow completion marker
        "Analysis/.workflow_complete"

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
    params:
        time=config.get("qc_complete_time", "00:10:00")
    resources:
        mem_mb=config.get("qc_complete_memory", 1000)
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

# This rule ensures trimming is completed before alignment starts
rule trimming_complete:
    input:
        trimmed=expand("Analysis/Trimmed/{sample}/{sample}_R1_trimmed.fastq.gz",
                      sample=SAMPLES.keys()),
        multiqc="Analysis/QC/Trimming/MultiQC/multiqc_report.html",
        multiqc_data="Analysis/QC/Trimming/MultiQC/multiqc_report_data"
    output:
        touch("Analysis/Trimmed/.trimming_complete")
    log:
        "logs/workflow/trimming_complete.log"
    params:
        time=config.get("trimming_complete_time", "00:10:00")
    resources:
        mem_mb=config.get("trimming_complete_memory", 1000)
    shell:
        """
        # Create logs directory if it doesn't exist
        mkdir -p $(dirname {log})
        
        # Log trimming completion
        echo "Trimming completed successfully at $(date)" | tee {log}
        echo "Trimmed samples: $(find Analysis/Trimmed -name '*_R1_trimmed.fastq.gz' | wc -l)" | tee -a {log}
        echo "MultiQC report: {input.multiqc}" | tee -a {log}
        
        # Verify that all required files exist
        if [ ! -f "{input.multiqc}" ]; then
            echo "ERROR: MultiQC report not found" | tee -a {log}
            exit 1
        fi
        
        # Verify that the multiqc_data directory exists
        if [ ! -d "{input.multiqc_data}" ]; then
            echo "WARNING: MultiQC data directory not found, creating placeholder" | tee -a {log}
            mkdir -p {input.multiqc_data}
            touch {input.multiqc_data}/.placeholder
        fi
        
        # Create a marker file to indicate trimming completion
        echo "Creating trimming completion marker file" | tee -a {log}
        """

# Rule to generate trimmed_samples.csv file
rule generate_trimmed_samples_csv:
    input:
        trimming_complete="Analysis/Trimmed/.trimming_complete"
    output:
        csv="resources/metadata/trimmed_samples.csv"
    log:
        "logs/workflow/generate_trimmed_samples.log"
    params:
        time=config.get("generate_trimmed_csv_time", "00:10:00")
    resources:
        mem_mb=config.get("generate_trimmed_csv_memory", 1000)
    shell:
        """
        # Create output directory if it doesn't exist
        mkdir -p $(dirname {output.csv})
        
        # Create header for the CSV file
        echo "sample,R1,R2" > {output.csv}
        
        # Debug information
        echo "===== GENERATING TRIMMED SAMPLES CSV =====" > {log} 2>&1
        echo "Date: $(date)" >> {log} 2>&1
        echo "Working directory: $(pwd)" >> {log} 2>&1
        echo "Output CSV: {output.csv}" >> {log} 2>&1
        
        # Find all trimmed R1 files and generate the CSV entries
        echo "Searching for trimmed files in Analysis/Trimmed/..." >> {log} 2>&1
        find Analysis/Trimmed -name "*_R1_trimmed.fastq.gz" | sort >> {log} 2>&1
        
        # Process each R1 file
        for r1_file in $(find Analysis/Trimmed -name "*_R1_trimmed.fastq.gz" | sort); do
            # Extract sample name from the path
            sample=$(basename $(dirname $r1_file))
            echo "Processing sample: $sample" >> {log} 2>&1
            
            # Get the corresponding R2 file
            r2_file=$(echo $r1_file | sed 's/_R1_/_R2_/')
            
            # Verify R2 file exists
            if [ ! -f "$r2_file" ]; then
                echo "ERROR: R2 file not found for $r1_file" | tee -a {log}
                exit 1
            fi
            
            # Verify files are not empty
            if [ ! -s "$r1_file" ]; then
                echo "ERROR: R1 file is empty: $r1_file" | tee -a {log}
                exit 1
            fi
            
            if [ ! -s "$r2_file" ]; then
                echo "ERROR: R2 file is empty: $r2_file" | tee -a {log}
                exit 1
            fi
            
            # Add entry to CSV file
            echo "$sample,$r1_file,$r2_file" >> {output.csv}
            echo "Added sample $sample to trimmed samples CSV" >> {log} 2>&1
            echo "  R1: $r1_file ($(ls -lh $r1_file | awk '{{print $5}}'))" >> {log} 2>&1
            echo "  R2: $r2_file ($(ls -lh $r2_file | awk '{{print $5}}'))" >> {log} 2>&1
        done
        
        # Verify the CSV file was created and has entries
        if [ ! -s "{output.csv}" ]; then
            echo "ERROR: CSV file is empty or was not created properly" | tee -a {log}
            exit 1
        fi
        
        # Count the number of samples in the CSV (excluding header)
        sample_count=$(wc -l < {output.csv})
        sample_count=$((sample_count - 1))  # Subtract 1 for the header
        
        if [ "$sample_count" -eq 0 ]; then
            echo "ERROR: No samples found in the CSV file" | tee -a {log}
            exit 1
        fi
        
        # Log completion
        echo "Generated trimmed samples CSV with $sample_count samples" | tee -a {log}
        echo "CSV file: {output.csv}" | tee -a {log}
        echo "CSV contents:" >> {log} 2>&1
        cat {output.csv} >> {log} 2>&1
        echo "===== CSV GENERATION COMPLETED =====" >> {log} 2>&1
        """

# Modify the original fastp rule to require QC completion
# This ensures proper stage dependencies
ruleorder: fastp_with_dependency > fastp

def get_fastp_params(wildcards, params_file="Analysis/QC/Trimming/trimming_params.json"):
    """
    Get fastp parameters for a sample directly in Python
    
    This function reads the trimming_params.json file generated by the generate_trimming_params rule
    in workflow/qc_params.snakefile. The JSON file contains sample-specific parameters determined
    by analyzing FastQC results.
    
    If the JSON file doesn't exist or doesn't contain parameters for the sample,
    default parameters are used.
    """
    import json
    import re
    import os
    
    # Default parameters
    defaults = {
        "leading_quality": 3,
        "trailing_quality": 3,
        "min_length": 36,
        "sliding_window": "4:15"
    }
    
    # Define large samples that need special handling
    large_samples = ["Scaber_SRR28516486", "Scaber_SRR28516488", "Scaber_SRR28516489"]
    
    # Load parameters file
    sample_params = {}
    try:
        if os.path.exists(params_file):
            with open(params_file) as f:
                all_params = json.load(f)
                
            # First try direct sample name match
            if wildcards.sample in all_params:
                sample_params = all_params[wildcards.sample]
                print(f"Found parameters for sample: {wildcards.sample}")
            else:
                # Try to extract SRR ID from the sample name and use that
                srr_match = re.search(r'(SRR\d+)', wildcards.sample)
                if srr_match and srr_match.group(1) in all_params:
                    srr_id = srr_match.group(1)
                    sample_params = all_params[srr_id]
                    print(f"Found parameters via SRR ID: {srr_id}")
                else:
                    print(f"No parameters found for {wildcards.sample}. Using defaults.")
                    sample_params = defaults
        else:
            print(f"Warning: Parameters file {params_file} not found, using defaults")
            sample_params = defaults
    except Exception as e:
        print(f"Error loading parameters: {str(e)}, using defaults")
        sample_params = defaults
    
    # Disable deduplication for large samples to prevent memory issues
    if wildcards.sample in large_samples and 'deduplication' in sample_params:
        print(f"Disabling deduplication for large sample: {wildcards.sample}")
        sample_params['deduplication'] = False
    
    # Build fastp command options
    cmd_parts = []
    
    # Quality filtering parameters
    cmd_parts.append(f"--qualified_quality_phred {sample_params.get('trailing_quality', 3)}")
    cmd_parts.append("--unqualified_percent_limit 40")
    
    # Sliding window parameters
    sliding_window = sample_params.get('sliding_window', '4:15')
    window_size, window_quality = sliding_window.split(':')
    cmd_parts.append("--cut_front")
    cmd_parts.append(f"--cut_front_window_size {window_size}")
    cmd_parts.append(f"--cut_front_mean_quality {window_quality}")
    cmd_parts.append("--cut_tail")
    cmd_parts.append(f"--cut_tail_window_size {window_size}")
    cmd_parts.append(f"--cut_tail_mean_quality {window_quality}")
    
    # Length filtering
    cmd_parts.append(f"--length_required {sample_params.get('min_length', 36)}")
    
    # Add optional flags
    if sample_params.get('deduplication', False):
        # For samples with deduplication, use a smaller hash size to reduce memory usage
        cmd_parts.append("--dedup")
        cmd_parts.append("--dedup_hash_size 24")  # Reduce hash size from default 32
    
    if sample_params.get('trim_poly_g', False):
        cmd_parts.append("--trim_poly_g")
    
    if sample_params.get('low_complexity_filter', False):
        cmd_parts.append("--low_complexity_filter")
    
    # Return the complete command string
    return " ".join(cmd_parts)

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
        fastp_opts=get_fastp_params,
        time=config.get("trimming_time", "02:00:00")
    threads: config.get("trimming_threads", 4)
    resources:
        mem_mb=config.get("trimming_memory", 32000)  # Increased from 8GB to 32GB
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

# Modify the star_align rule to require trimming completion
rule star_align_with_dependency:
    input:
        # We'll keep these as dependencies but won't use them directly
        # Instead, we'll read the actual paths from the CSV file
        r1="Analysis/Trimmed/{sample}/{sample}_R1_trimmed.fastq.gz",
        r2="Analysis/Trimmed/{sample}/{sample}_R2_trimmed.fastq.gz",
        index=config.get("star_index", "resources/genome/star_index"),
        trimming_complete="Analysis/Trimmed/.trimming_complete",
        trimmed_samples_csv="resources/metadata/trimmed_samples.csv"
    output:
        bam="Analysis/Alignment/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
        transcriptome_bam="Analysis/Alignment/STAR/{sample}/{sample}.Aligned.toTranscriptome.out.bam",
        counts="Analysis/Alignment/STAR/{sample}/{sample}.ReadsPerGene.out.tab",
        sj="Analysis/Alignment/STAR/{sample}/{sample}.SJ.out.tab",
        log_final="Analysis/Alignment/STAR/{sample}/{sample}.Log.final.out",
        log="Analysis/Alignment/STAR/{sample}/{sample}.Log.out",
        log_progress="Analysis/Alignment/STAR/{sample}/{sample}.Log.progress.out"
    log:
        "logs/star/{sample}.log"
    params:
        # STAR alignment parameters
        gtf=config.get("gtf_file", "resources/genome/Homo_sapiens.GRCh38.90.gtf"),
        prefix="Analysis/Alignment/STAR/{sample}/{sample}.",
        # STAR specific parameters
        outFilterMismatchNmax=config.get("star_mismatch", 10),
        outFilterMultimapNmax=config.get("star_multimap", 20),
        outFilterScoreMinOverLread=config.get("star_score_min", 0.66),
        outFilterMatchNminOverLread=config.get("star_match_min", 0.66),
        alignSJDBoverhangMin=config.get("star_sjdb_overhang_min", 3),
        alignIntronMax=config.get("star_intron_max", 500000),
        # SLURM parameters
        time=config.get("star_align_time", "04:00:00")
    threads: config.get("star_threads", 8)
    resources:
        mem_mb=config.get("star_memory", 32000)  # 32GB
    shell:
        """
        # Create output directory
        mkdir -p $(dirname {output.bam})
        
        # Log input files and verify they exist
        echo "===== STAR ALIGNMENT STARTING =====" > {log} 2>&1
        echo "Processing sample: {wildcards.sample}" >> {log} 2>&1
        echo "Date: $(date)" >> {log} 2>&1
        echo "Hostname: $(hostname)" >> {log} 2>&1
        echo "Working directory: $(pwd)" >> {log} 2>&1
        
        # Get the actual trimmed file paths from the CSV
        echo "Reading trimmed file paths from CSV..." >> {log} 2>&1
        if [ ! -f "{input.trimmed_samples_csv}" ]; then
            echo "ERROR: Trimmed samples CSV not found: {input.trimmed_samples_csv}" >> {log} 2>&1
            exit 1
        fi
        
        # Extract the R1 and R2 paths for this sample from the CSV
        CSV_LINE=$(grep "^{wildcards.sample}," {input.trimmed_samples_csv})
        if [ -z "$CSV_LINE" ]; then
            echo "ERROR: Sample {wildcards.sample} not found in {input.trimmed_samples_csv}" >> {log} 2>&1
            exit 1
        fi
        
        # Parse the CSV line to get R1 and R2 paths
        R1_PATH=$(echo "$CSV_LINE" | cut -d',' -f2)
        R2_PATH=$(echo "$CSV_LINE" | cut -d',' -f3)
        
        echo "Using trimmed files from CSV:" >> {log} 2>&1
        echo "  R1: $R1_PATH" >> {log} 2>&1
        echo "  R2: $R2_PATH" >> {log} 2>&1
        
        # Verify the files exist
        if [ ! -f "$R1_PATH" ]; then
            echo "ERROR: R1 file from CSV not found: $R1_PATH" >> {log} 2>&1
            exit 1
        fi
        
        if [ ! -f "$R2_PATH" ]; then
            echo "ERROR: R2 file from CSV not found: $R2_PATH" >> {log} 2>&1
            exit 1
        fi
        
        # Verify STAR index
        echo "Checking STAR index..." >> {log} 2>&1
        if [ ! -d "{input.index}" ]; then
            echo "ERROR: STAR index directory not found: {input.index}" >> {log} 2>&1
            exit 1
        else
            echo "STAR index: {input.index}" >> {log} 2>&1
            echo "Index files:" >> {log} 2>&1
            ls -lh {input.index} >> {log} 2>&1
        fi
        
        # Verify GTF file
        echo "Checking GTF file..." >> {log} 2>&1
        if [ ! -f "{params.gtf}" ]; then
            echo "ERROR: GTF file not found: {params.gtf}" >> {log} 2>&1
            exit 1
        else
            echo "GTF file: {params.gtf} ($(ls -lh {params.gtf} | awk '{{print $5}}'))" >> {log} 2>&1
        fi
        
        # Check available memory
        echo "Available memory:" >> {log} 2>&1
        free -h >> {log} 2>&1
        
        # Check STAR version
        echo "STAR version:" >> {log} 2>&1
        STAR --version >> {log} 2>&1
        
        echo "===== STARTING STAR ALIGNMENT =====" >> {log} 2>&1
        
        # Run STAR alignment with paths from CSV
        STAR --runThreadN {threads} \\
            --genomeDir {input.index} \\
            --readFilesIn $R1_PATH $R2_PATH \\
            --readFilesCommand zcat \\
            --outFileNamePrefix {params.prefix} \\
            --outSAMtype BAM SortedByCoordinate \\
            --outSAMunmapped Within \\
            --outSAMattributes Standard \\
            --outFilterMismatchNmax {params.outFilterMismatchNmax} \\
            --outFilterMultimapNmax {params.outFilterMultimapNmax} \\
            --outFilterScoreMinOverLread {params.outFilterScoreMinOverLread} \\
            --outFilterMatchNminOverLread {params.outFilterMatchNminOverLread} \\
            --alignSJDBoverhangMin {params.alignSJDBoverhangMin} \\
            --alignIntronMax {params.alignIntronMax} \\
            --quantMode TranscriptomeSAM GeneCounts \\
            --sjdbGTFfile {params.gtf} \\
            >> {log} 2>&1
        
        # Check STAR exit status
        STAR_EXIT=$?
        if [ $STAR_EXIT -ne 0 ]; then
            echo "ERROR: STAR alignment failed with exit code $STAR_EXIT" >> {log} 2>&1
            exit $STAR_EXIT
        fi
        
        echo "===== VERIFYING OUTPUTS =====" >> {log} 2>&1
        
        # Verify output files exist and have non-zero size
        for output_file in "{output.bam}" "{output.transcriptome_bam}" "{output.counts}" "{output.sj}" "{output.log_final}" "{output.log}" "{output.log_progress}"; do
            if [ ! -f "$output_file" ]; then
                echo "ERROR: Output file not created: $output_file" >> {log} 2>&1
                exit 1
            elif [ ! -s "$output_file" ]; then
                echo "ERROR: Output file is empty: $output_file" >> {log} 2>&1
                exit 1
            else
                echo "Output file created: $output_file ($(ls -lh $output_file | awk '{{print $5}}'))" >> {log} 2>&1
            fi
        done
        
        # Check BAM file with samtools
        echo "Checking BAM file with samtools..." >> {log} 2>&1
        samtools quickcheck {output.bam}
        if [ $? -ne 0 ]; then
            echo "ERROR: BAM file failed samtools quickcheck: {output.bam}" >> {log} 2>&1
            exit 1
        fi
        
        echo "===== STAR ALIGNMENT COMPLETED SUCCESSFULLY =====" >> {log} 2>&1
        echo "Date: $(date)" >> {log} 2>&1
        """

# Set rule order for alignment with dependency
ruleorder: star_align_with_dependency > star_align
ruleorder: main_alignment_complete > alignment_complete

# This rule ensures alignment is completed
rule main_alignment_complete:
    input:
        bams=expand("Analysis/Alignment/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
               sample=SAMPLES.keys()),
        bai=expand("Analysis/Alignment/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai",
               sample=SAMPLES.keys()),
        multiqc="Analysis/Alignment/MultiQC/multiqc_report.html"
    output:
        touch("Analysis/Alignment/.main_alignment_complete")
    log:
        "logs/workflow/main_alignment_complete.log"
    params:
        time=config.get("alignment_complete_time", "00:10:00")
    resources:
        mem_mb=config.get("alignment_complete_memory", 1000)
    shell:
        """
        # Create logs directory if it doesn't exist
        mkdir -p $(dirname {log})
        
        # Log alignment completion
        echo "Alignment completed successfully at $(date)" | tee {log}
        echo "Aligned samples: $(find Analysis/Alignment/STAR -name '*.Aligned.sortedByCoord.out.bam' | wc -l)" | tee -a {log}
        echo "MultiQC report: {input.multiqc}" | tee -a {log}
        
        # Verify that all required files exist
        if [ ! -f "{input.multiqc}" ]; then
            echo "ERROR: MultiQC report not found" | tee -a {log}
            exit 1
        fi
        
        # Create a marker file to indicate alignment completion
        echo "Creating alignment completion marker file" | tee -a {log}
        
        # Also create the standard alignment complete marker for compatibility
        touch Analysis/Alignment/.alignment_complete
        """

# Rules for running specific workflow stages independently
# These allow users to run just a portion of the workflow

# Run just the parameters generation step
rule params_generation:
    input:
        "Analysis/QC/Trimming/trimming_params.json"
    params:
        time=config.get("params_generation_time", "00:05:00")
    resources:
        mem_mb=config.get("params_generation_memory", 1000)
    shell:
        """
        echo "Parameter generation completed successfully"
        """

# Run just the QC stage
rule qc_stage:
    input:
        "Analysis/QC/.qc_complete"
    params:
        time=config.get("qc_stage_time", "00:05:00")
    resources:
        mem_mb=config.get("qc_stage_memory", 1000)
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
        "Analysis/QC/Trimming/MultiQC/multiqc_report.html",
        "Analysis/Trimmed/.trimming_complete",
        "resources/metadata/trimmed_samples.csv"
    params:
        time=config.get("trimming_stage_time", "00:05:00")
    resources:
        mem_mb=config.get("trimming_stage_memory", 1000)
    shell:
        """
        echo "Trimming stage completed successfully"
        echo "Trimmed samples CSV generated at resources/metadata/trimmed_samples.csv"
        """

# Run just the alignment stage
rule alignment_stage:
    input:
        "Analysis/Alignment/.main_alignment_complete"
    params:
        time=config.get("alignment_stage_time", "00:05:00")
    resources:
        mem_mb=config.get("alignment_stage_memory", 1000)
    shell:
        """
        echo "Alignment stage completed successfully"
        """

# Define workflow stages for better control over execution order
rule workflow_stages:
    input:
        qc="Analysis/QC/.qc_complete",
        trimming="Analysis/Trimmed/.trimming_complete",
        alignment="Analysis/Alignment/.main_alignment_complete"
    output:
        touch("Analysis/.workflow_complete")
    params:
        time=config.get("workflow_stages_time", "00:10:00")
    resources:
        mem_mb=config.get("workflow_stages_memory", 1000)
    shell:
        """
        echo "Workflow completed successfully at $(date)"
        echo "QC stage: Completed"
        echo "Trimming stage: Completed"
        echo "Alignment stage: Completed"
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
    print("Alignment BAMs: Analysis/Alignment/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam")
    print("Alignment counts: Analysis/Alignment/STAR/{sample}/{sample}.ReadsPerGene.out.tab")
    print("Alignment MultiQC: Analysis/Alignment/MultiQC/multiqc_report.html")
    print("\nTo visualize the workflow graph:")
    print("  snakemake --dag | dot -Tsvg > workflow.svg")

onerror:
    print("\nWorkflow failed!")
    print("===============")
    print("Check log files in the logs/ directory for error details.")
    print("Run with --notemp flag to keep temporary files for debugging.")
    print("Use -n (dry run) and -p (print shell commands) for troubleshooting.")