from snakemake.utils import min_version
min_version("6.0")

# Define SAMPLES as an empty dict if not already defined (for syntax checking)
if 'SAMPLES' not in globals():
    SAMPLES = {}
    config = {}

# No configfile import - config is loaded in main Snakefile

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
        sjdb_overhang=config.get("sjdb_overhang", 100)  # Read length - 1
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

# Rule to index BAM files
rule index_bam:
    input:
        bam="Analysis/Alignment/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        bai="Analysis/Alignment/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai"
    log:
        "logs/star/{sample}.index.log"
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
        outdir="Analysis/Alignment/MultiQC"
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
        
        # Create data directory if it doesn't exist
        if [ ! -d "{output.data_dir}" ]; then
            mkdir -p {output.data_dir}
            touch {output.data_dir}/.placeholder
        fi
        """