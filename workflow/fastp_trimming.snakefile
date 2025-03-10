from snakemake.utils import min_version
min_version("6.0")

# No configfile import - config is loaded in main Snakefile

# Rule to run trimming for all samples
rule trim_all:
    input:
        trimmed=expand("Analysis/Trimmed/{sample}/{sample}_R{read}_trimmed.fastq.gz",
                     sample=SAMPLES.keys(),
                     read=["1", "2"]),
        report="Analysis/QC/Trimming/MultiQC/multiqc_report.html"

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