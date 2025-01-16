import os

# Create the output directory
os.makedirs("Analysis/QC/MultiQC", exist_ok=True)
os.makedirs("logs/multiqc", exist_ok=True)

# Simple rule to run MultiQC
rule all:
    input:
        "Analysis/QC/MultiQC/multiqc_report.html",
        "Analysis/QC/MultiQC/multiqc_report.pdf"

rule multiqc:
    input:
        # Use the directory containing all FastQC outputs
        fastqc_dir="Analysis/QC/FastQC"
    output:
        html="Analysis/QC/MultiQC/multiqc_report.html",
        pdf="Analysis/QC/MultiQC/multiqc_report.pdf"
    log:
        "logs/multiqc/multiqc.log"
    params:
        outdir="Analysis/QC/MultiQC"
    resources:
        mem_mb=4000,
        time="00:30:00"
    shell:
        """
        # Run MultiQC
        multiqc --force \
            --outdir {params.outdir} \
            {input.fastqc_dir} \
            2> {log}
            
        # Convert to PDF
        wkhtmltopdf {output.html} {output.pdf}
        """