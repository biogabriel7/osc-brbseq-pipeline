import os
import csv
from snakemake.utils import min_version
min_version("6.0")

# Load configuration
configfile: "resources/config/params.yaml"

# Create necessary directories
for dir_path in [
    "logs/star", 
    "Analysis/Alignment/STAR", 
    "Analysis/Alignment/MultiQC"
]:
    os.makedirs(dir_path, exist_ok=True)

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
    except FileNotFoundError:
        print(f"Error: Metasheet not found at {config['metasheet_path']}")
        raise
    except KeyError as e:
        print(f"Error: Missing required column in metasheet: {e}")
        raise
    return samples

# Store samples from metasheet in a dictionary
SAMPLES = get_samples_from_metasheet()

# Rule to run alignment for all samples
rule align_all:
    input:
        expand("Analysis/Alignment/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai",
               sample=SAMPLES.keys()),
        "Analysis/Alignment/MultiQC/multiqc_report.html"

# Rule to create STAR index (run this once before the main workflow)
rule create_star_index:
    input:
        fasta=config.get("genome_fasta", "resources/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa"),
        gtf=config.get("gtf_file", "resources/genome/Homo_sapiens.GRCh38.90.gtf")
    output:
        directory(config.get("star_index", "resources/genome/star_index"))
    log:
        "logs/star/star_index.log"
    params:
        sjdb_overhang=config.get("sjdb_overhang", 100),  # Read length - 1
        time=config.get("star_index_time", "04:00:00")
    threads: config.get("star_index_threads", 16)
    resources:
        mem_mb=config.get("star_index_memory", 40000)  # 40GB
    shell:
        """
        mkdir -p {output}
        
        STAR --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeDir {output} \
            --genomeFastaFiles {input.fasta} \
            --sjdbGTFfile {input.gtf} \
            --sjdbOverhang {params.sjdb_overhang} \
            > {log} 2>&1
        """

# Rule to run STAR alignment with dependency on trimming
rule star_align:
    input:
        r1="Analysis/Trimmed/{sample}/{sample}_R1_trimmed.fastq.gz",
        r2="Analysis/Trimmed/{sample}/{sample}_R2_trimmed.fastq.gz",
        index=config.get("star_index", "resources/genome/star_index")
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
        
        # Run STAR alignment
        STAR --runThreadN {threads} \
            --genomeDir {input.index} \
            --readFilesIn {input.r1} {input.r2} \
            --readFilesCommand zcat \
            --outFileNamePrefix {params.prefix} \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMunmapped Within \
            --outSAMattributes Standard \
            --outFilterMismatchNmax {params.outFilterMismatchNmax} \
            --outFilterMultimapNmax {params.outFilterMultimapNmax} \
            --outFilterScoreMinOverLread {params.outFilterScoreMinOverLread} \
            --outFilterMatchNminOverLread {params.outFilterMatchNminOverLread} \
            --alignSJDBoverhangMin {params.alignSJDBoverhangMin} \
            --alignIntronMax {params.alignIntronMax} \
            --quantMode TranscriptomeSAM GeneCounts \
            --sjdbGTFfile {params.gtf} \
            > {log} 2>&1
        
        # Verify output files exist
        if [ ! -f "{output.bam}" ]; then
            echo "ERROR: BAM file was not created" >> {log}
            exit 1
        fi
        """

# Rule to index BAM files
rule index_bam:
    input:
        bam="Analysis/Alignment/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        bai="Analysis/Alignment/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai"
    log:
        "logs/star/{sample}.index.log"
    params:
        time=config.get("bam_index_time", "01:00:00")
    threads: 1
    resources:
        mem_mb=config.get("bam_index_memory", 4000)
    shell:
        """
        samtools index {input.bam} {output.bai} 2> {log}
        """

# Rule to extract alignment metrics
rule alignment_metrics:
    input:
        bam="Analysis/Alignment/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
        bai="Analysis/Alignment/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai"
    output:
        flagstat="Analysis/Alignment/STAR/{sample}/{sample}.flagstat.txt",
        idxstats="Analysis/Alignment/STAR/{sample}/{sample}.idxstats.txt",
        stats="Analysis/Alignment/STAR/{sample}/{sample}.stats.txt"
    log:
        "logs/star/{sample}.metrics.log"
    params:
        time=config.get("metrics_time", "01:00:00")
    threads: 1
    resources:
        mem_mb=config.get("metrics_memory", 4000)
    shell:
        """
        # Generate alignment metrics
        samtools flagstat {input.bam} > {output.flagstat} 2> {log}
        samtools idxstats {input.bam} > {output.idxstats} 2>> {log}
        samtools stats {input.bam} > {output.stats} 2>> {log}
        """

# Rule to run MultiQC on alignment results
rule multiqc_alignment:
    input:
        star_logs=expand("Analysis/Alignment/STAR/{sample}/{sample}.Log.final.out", sample=SAMPLES.keys()),
        flagstats=expand("Analysis/Alignment/STAR/{sample}/{sample}.flagstat.txt", sample=SAMPLES.keys()),
        idxstats=expand("Analysis/Alignment/STAR/{sample}/{sample}.idxstats.txt", sample=SAMPLES.keys()),
        stats=expand("Analysis/Alignment/STAR/{sample}/{sample}.stats.txt", sample=SAMPLES.keys())
    output:
        html="Analysis/Alignment/MultiQC/multiqc_report.html",
        data_dir=directory("Analysis/Alignment/MultiQC/multiqc_data")
    log:
        "logs/star/multiqc.log"
    params:
        outdir="Analysis/Alignment/MultiQC",
        time=config.get("multiqc_alignment_time", "01:00:00")
    resources:
        mem_mb=config.get("multiqc_alignment_memory", 4000)
    shell:
        """
        # Create output directory
        mkdir -p {params.outdir}
        
        # Run MultiQC
        multiqc --force \
            --outdir {params.outdir} \
            Analysis/Alignment/STAR/ \
            > {log} 2>&1
            
        # Verify output
        if [ ! -f "{output.html}" ]; then
            echo "ERROR: MultiQC report was not created" >> {log}
            exit 1
        fi
        """

# Rule to mark alignment completion
rule alignment_complete:
    input:
        bams=expand("Analysis/Alignment/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai", sample=SAMPLES.keys()),
        multiqc="Analysis/Alignment/MultiQC/multiqc_report.html"
    output:
        touch("Analysis/Alignment/.alignment_complete")
    log:
        "logs/workflow/alignment_complete.log"
    params:
        time=config.get("alignment_complete_time", "00:10:00")
    resources:
        mem_mb=config.get("alignment_complete_memory", 1000)
    shell:
        """
        # Log alignment completion
        echo "Alignment completed successfully at $(date)" | tee {log}
        echo "Aligned samples: $(find Analysis/Alignment/STAR -name '*.Aligned.sortedByCoord.out.bam' | wc -l)" | tee -a {log}
        echo "MultiQC report: {input.multiqc}" | tee -a {log}
        
        # Create a marker file to indicate alignment completion
        echo "Creating alignment completion marker file" | tee -a {log}
        """
