import csv
import yaml
from snakemake.utils import min_version
min_version("6.0")

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
                    "R2": row["R2"]
                }
                print(f"Sample: {sample_name}, R1: {samples[sample_name]['R1']}, R2: {samples[sample_name]['R2']}")
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
for dir in ["logs/fastqc", "Analysis/QC/FastQC", "Analysis/QC/MultiQC"]:
    os.makedirs(dir, exist_ok=True)

# Rule to generate all QC reports including MultiQC
rule all:
    input:
        expand("Analysis/QC/FastQC/{sample}/{sample}_{read}_fastqc.html",
               sample=SAMPLES.keys(), read=["R1", "R2"]),
        "Analysis/QC/MultiQC/multiqc_report.html",
        "Analysis/QC/MultiQC/multiqc_report.pdf"

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
    threads: config["fastqc_threads"] if "fastqc_threads" in config else 2
    resources:
        mem_mb=config["fastqc_memory"] if "fastqc_memory" in config else 4000,
        time=config["fastqc_time"] if "fastqc_time" in config else "01:00:00"
    shell:
        """
        mkdir -p {params.outdir}
        module load fastqc
        
        fastqc --threads {threads} \
            --outdir {params.outdir} \
            {input.r1} {input.r2} \
            2> {log}
        """

# Rule to run MultiQC to aggregate all FastQC reports
rule multiqc:
    input:
        expand("Analysis/QC/FastQC/{sample}/{sample}_{read}_fastqc.zip",
               sample=SAMPLES.keys(), read=["R1", "R2"])
    output:
        html="Analysis/QC/MultiQC/multiqc_report.html",
        pdf="Analysis/QC/MultiQC/multiqc_report.pdf"
    log:
        "logs/fastqc/multiqc.log"
    params:
        outdir="Analysis/QC/MultiQC"
    resources:
        mem_mb=config["multiqc_memory"] if "multiqc_memory" in config else 4000,
        time=config["multiqc_time"] if "multiqc_time" in config else "00:30:00"
    shell:
        """
        module load python/3.8-conda
        source activate multiqc_env
        
        multiqc --force \
            --outdir {params.outdir} \
            Analysis/QC/FastQC/ \
            2> {log}
            
        # Convert HTML to PDF using wkhtmltopdf
        module load wkhtmltopdf
        wkhtmltopdf {output.html} {output.pdf}
        """