from snakemake.utils import min_version
import os
import json
import re

# Ensure minimum Snakemake version for compatibility
min_version("6.0")

# Import centralized utility functions
from workflow.common.utils import get_samples_from_metasheet, create_workflow_dirs

# Create essential directories
create_workflow_dirs()

# Load configuration file (using relative path for portability)
configfile: "resources/config/params.yaml"

# Check if we're in test mode
TEST_MODE = config.get("test_mode", False)
if TEST_MODE:
    print("RUNNING IN TEST MODE WITH REDUCED DATASET")

# Get samples information
SAMPLES = get_samples_from_metasheet(config, test_mode=TEST_MODE)

# Include rules from other Snakefiles
# This brings in all the rules defined in the other files
include: "workflow/qc_params.snakefile"      # QC and parameter generation rules
include: "workflow/fastp_trimming.snakefile" # Trimming rules
include: "workflow/star_alignment.snakefile" # Alignment rules

# Define the workflow order
# Updated to resolve ambiguity between fastp and fastp_trim
ruleorder: fastp > fastp_trim > fastqc > multiqc > generate_trimming_params > multiqc_fastp > star_align > index_bam > alignment_metrics > multiqc_alignment

# Define the complete workflow with all expected outputs
rule all:
    input:
        # QC stage outputs
        expand("Analysis/QC/FastQC/{sample}/{sample}_R{read}_fastqc.html",
               sample=SAMPLES.keys(),
               read=["1", "2"]),
        "Analysis/QC/MultiQC/multiqc_report.html",
        "Analysis/QC/Trimming/trimming_params.json",
        # QC completion marker
        "Analysis/QC/.qc_complete",
        # Trimming stage outputs
        expand("Analysis/Trimmed/{sample}/{sample}_R{read}_trimmed.fastq.gz",
               sample=SAMPLES.keys(),
               read=["1", "2"]),
        "Analysis/QC/Trimming/MultiQC/multiqc_report.html",
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
    shell:
        """
        echo "QC completed successfully at $(date)" > {log}
        """

# This rule ensures trimming is completed before alignment starts
rule trimming_complete:
    input:
        trimmed=expand("Analysis/Trimmed/{sample}/{sample}_R{read}_trimmed.fastq.gz",
                      sample=SAMPLES.keys(),
                      read=["1", "2"]),
        multiqc="Analysis/QC/Trimming/MultiQC/multiqc_report.html"
    output:
        touch("Analysis/Trimmed/.trimming_complete")
    log:
        "logs/workflow/trimming_complete.log"
    shell:
        """
        echo "Trimming completed successfully at $(date)" > {log}
        """

# Rule to generate trimmed_samples.csv file
rule generate_trimmed_samples_csv:
    input:
        trimming_complete="Analysis/Trimmed/.trimming_complete"
    output:
        csv="resources/metadata/trimmed_samples.csv"
    log:
        "logs/workflow/generate_trimmed_samples.log"
    shell:
        """
        # Create header for the CSV file
        echo "sample,R1,R2" > {output.csv}
        
        # Process each R1 file
        for r1_file in $(find Analysis/Trimmed -name "*_R1_trimmed.fastq.gz" | sort); do
            # Extract sample name from the path
            sample=$(basename $(dirname $r1_file))
            
            # Get the corresponding R2 file
            r2_file=$(echo $r1_file | sed 's/_R1_/_R2_/')
            
            # Add entry to CSV file
            echo "$sample,$r1_file,$r2_file" >> {output.csv}
            echo "Added sample $sample to trimmed samples CSV" >> {log}
        done
        
        echo "Generated trimmed samples CSV with $(grep -c ',' {output.csv} | awk '{{print $1-1}}') samples" | tee -a {log}
        """

# Get fastp parameters function for the fastp rule
def get_fastp_params(wildcards, params_file="Analysis/QC/Trimming/trimming_params.json"):
    """
    Get fastp parameters for a sample
    
    This function reads the trimming_params.json file generated by the generate_trimming_params rule
    and returns sample-specific parameters.
    """
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
    
    if sample_params.get('trim_poly_g', False):
        cmd_parts.append("--trim_poly_g")
    
    if sample_params.get('low_complexity_filter', False):
        cmd_parts.append("--low_complexity_filter")
    
    # Return the complete command string
    return " ".join(cmd_parts)

# Primary fastp rule
rule fastp:
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
        fastp_opts=get_fastp_params
    threads: config.get("trimming_threads", 4)
    resources:
        mem_mb=config.get("trimming_memory", 32000)
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
        set -o pipefail  # Ensure pipeline errors are caught
        fastp \\
            --in1 {input.r1} \\
            --in2 {input.r2} \\
            --out1 {output.r1} \\
            --out2 {output.r2} \\
            --html {output.html} \\
            --json {output.json} \\
            --thread {threads} \\
            --detect_adapter_for_pe \\
            {params.fastp_opts} \\
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

# STAR alignment rule with dependency on trimming completion
rule star_align:
    input:
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
        alignIntronMax=config.get("star_intron_max", 500000)
    threads: config.get("star_threads", 8)
    resources:
        mem_mb=config.get("star_memory", 32000)
    shell:
        """
        # Create output directory
        mkdir -p $(dirname {output.bam})
        
        # Log input files and verify they exist
        echo "===== STAR ALIGNMENT STARTING =====" > {log} 2>&1
        echo "Processing sample: {wildcards.sample}" >> {log} 2>&1
        echo "Date: $(date)" >> {log} 2>&1
        echo "Input reads:" >> {log} 2>&1
        echo "  R1: {input.r1}" >> {log} 2>&1
        echo "  R2: {input.r2}" >> {log} 2>&1
        echo "STAR index: {input.index}" >> {log} 2>&1
        echo "GTF file: {params.gtf}" >> {log} 2>&1
        echo "Threads: {threads}" >> {log} 2>&1
        
        # Verify input files exist
        if [ ! -f "{input.r1}" ] || [ ! -f "{input.r2}" ]; then
            echo "ERROR: Input fastq files not found" >> {log} 2>&1
            exit 1
        fi
        
        if [ ! -d "{input.index}" ]; then
            echo "ERROR: STAR index directory not found" >> {log} 2>&1
            exit 1
        fi
        
        if [ ! -f "{params.gtf}" ]; then
            echo "ERROR: GTF file not found" >> {log} 2>&1
            exit 1
        fi
        
        # Set temporary directory for STAR (helps with network file system issues)
        export TMPDIR=$(dirname {output.bam})/tmp
        mkdir -p $TMPDIR
        
        # Run STAR alignment 
        echo "Running STAR alignment..." >> {log} 2>&1
        set -o pipefail  # Ensure pipeline errors are caught
        
        STAR --runThreadN {threads} \\
            --genomeDir {input.index} \\
            --readFilesIn {input.r1} {input.r2} \\
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
            --limitBAMsortRAM $(({resources.mem_mb} * 1000000 * 3/4)) \\
            >> {log} 2>&1
        
        STAR_EXIT_CODE=$?
        if [ $STAR_EXIT_CODE -ne 0 ]; then
            echo "ERROR: STAR alignment failed with exit code $STAR_EXIT_CODE" >> {log} 2>&1
            exit $STAR_EXIT_CODE
        fi
        
        # Verify output files exist
        if [ ! -f "{output.bam}" ]; then
            echo "ERROR: BAM file was not created: {output.bam}" >> {log} 2>&1
            exit 1
        fi
        
        if [ ! -f "{output.counts}" ]; then
            echo "ERROR: Counts file was not created: {output.counts}" >> {log} 2>&1
            exit 1
        fi
        
        # Check BAM file with samtools
        samtools quickcheck {output.bam}
        if [ $? -ne 0 ]; then
            echo "ERROR: BAM file failed samtools quickcheck: {output.bam}" >> {log} 2>&1
            exit 1
        fi
        
        # Print alignment statistics
        echo "===== ALIGNMENT STATISTICS =====" >> {log} 2>&1
        grep "Number of input reads" {output.log_final} >> {log} 2>&1
        grep "Uniquely mapped reads" {output.log_final} >> {log} 2>&1
        grep "% of reads mapped to multiple loci" {output.log_final} >> {log} 2>&1
        grep "% of reads unmapped" {output.log_final} >> {log} 2>&1
        
        # Clean up temporary directory
        rm -rf $TMPDIR
        
        echo "===== STAR ALIGNMENT COMPLETED SUCCESSFULLY =====" >> {log} 2>&1
        echo "Date: $(date)" >> {log} 2>&1
        echo "Output BAM: {output.bam}" >> {log} 2>&1
        echo "Output counts: {output.counts}" >> {log} 2>&1
        """

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
    shell:
        """
        echo "Alignment completed successfully at $(date)" > {log}
        
        # Also create the standard alignment complete marker for compatibility
        touch Analysis/Alignment/.alignment_complete
        """

# Define workflow stages for better control over execution order
rule workflow_stages:
    input:
        qc="Analysis/QC/.qc_complete",
        trimming="Analysis/Trimmed/.trimming_complete",
        alignment="Analysis/Alignment/.main_alignment_complete"
    output:
        touch("Analysis/.workflow_complete")
    shell:
        """
        echo "Workflow completed successfully at $(date)"
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