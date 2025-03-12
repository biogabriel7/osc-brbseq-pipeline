from snakemake.utils import min_version
min_version("6.0")

# No configfile import - config is loaded in main Snakefile

# Define SAMPLES as an empty dict if not already defined (for syntax checking)
if 'SAMPLES' not in globals():
    SAMPLES = {}
    config = {}

# Rule to run trimming for all samples
rule trim_all:
    input:
        trimmed=expand("Analysis/Trimmed/{sample}/{sample}_R{read}_trimmed.fastq.gz",
                     sample=SAMPLES.keys(),
                     read=["1", "2"]),
        report="Analysis/QC/Trimming/MultiQC/multiqc_report.html"

# Rule to run fastp trimming for each sample
rule fastp_trim:
    input:
        r1=lambda wildcards: SAMPLES[wildcards.sample]["R1"],
        r2=lambda wildcards: SAMPLES[wildcards.sample]["R2"]
    output:
        r1="Analysis/Trimmed/{sample}/{sample}_R1_trimmed.fastq.gz",
        r2="Analysis/Trimmed/{sample}/{sample}_R2_trimmed.fastq.gz",
        json="Analysis/QC/Trimming/Reports/{sample}_fastp.json",
        html="Analysis/QC/Trimming/Reports/{sample}_fastp.html"
    log:
        "logs/trimming/{sample}_fastp.log"
    threads: config.get("fastp_threads", 4)
    resources:
        mem_mb=config.get("fastp_memory", 8000)
    shell:
        """
        # Create output directories
        mkdir -p $(dirname {output.r1})
        mkdir -p $(dirname {output.json})
        
        # Run fastp
        fastp \
            --in1 {input.r1} \
            --in2 {input.r2} \
            --out1 {output.r1} \
            --out2 {output.r2} \
            --json {output.json} \
            --html {output.html} \
            --thread {threads} \
            --detect_adapter_for_pe \
            --qualified_quality_phred 15 \
            --unqualified_percent_limit 40 \
            --n_base_limit 5 \
            --length_required 36 \
            --correction \
            --trim_front1 3 \
            --trim_front2 3 \
            > {log} 2>&1
        """

# Rule to run MultiQC on fastp reports
rule multiqc_fastp:
    input:
        fastp_reports=expand("Analysis/QC/Trimming/Reports/{sample}_fastp.json",
                          sample=SAMPLES.keys())
    output:
        html="Analysis/QC/Trimming/MultiQC/multiqc_report.html",
        data_dir=directory("Analysis/QC/Trimming/MultiQC/multiqc_report_data")
    log:
        "logs/trimming/multiqc_fastp.log"
    resources:
        mem_mb=config.get("multiqc_fastp_memory", 4000)
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
        
        # Remove any existing output directory to avoid conflicts
        rm -rf {output.data_dir}
        
        multiqc \
            --force \
            --outdir $(dirname {output.html}) \
            --filename $(basename {output.html}) \
            Analysis/QC/Trimming/Reports/ \
            >> {log} 2>&1
            
        # Verify outputs
        if [ ! -f "{output.html}" ]; then
            echo "ERROR: MultiQC report was not created" >> {log}
            exit 1
        fi
        
        # Create the data directory if it doesn't exist (sometimes MultiQC doesn't create it)
        if [ ! -d "{output.data_dir}" ]; then
            echo "Creating data directory manually" >> {log}
            mkdir -p {output.data_dir}
            touch {output.data_dir}/.placeholder
        fi
        
        echo "MultiQC completed successfully" >> {log}
        """