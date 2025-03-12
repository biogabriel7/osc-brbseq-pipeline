# RNA-seq Alignment

This document details the Alignment stage of the RNA-seq analysis workflow.

## Overview

The Alignment stage maps trimmed reads to a reference genome using STAR, generates alignment metrics, and creates a comprehensive report. This stage consists of four main steps:

1. **STAR Alignment**: Maps trimmed reads to the reference genome
2. **BAM Indexing**: Creates indices for the aligned BAM files
3. **Alignment Metrics**: Generates statistics about the alignments
4. **MultiQC Reporting**: Aggregates alignment metrics into a comprehensive report

## Workflow Steps

### 1. STAR Alignment

STAR aligns trimmed reads to the reference genome, using the automatically generated trimmed samples metadata file.

**Rule**: `star_align_with_dependency`

```python
rule star_align_with_dependency:
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
        alignIntronMax=config.get("star_intron_max", 500000),
        # SLURM parameters
        time=config.get("star_align_time", "04:00:00")
    threads: config.get("star_threads", 8)
    resources:
        mem_mb=config.get("star_memory", 32000)  # 32GB
```

Note that this rule requires the `trimmed_samples_csv` file as input, ensuring that the trimmed samples metadata is available before alignment begins.

### 2. BAM Indexing

This step creates indices for the aligned BAM files, which are required for many downstream analyses.

**Rule**: `index_bam`

```python
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
```

### 3. Alignment Metrics

This step generates various statistics about the alignments, which are useful for quality assessment.

**Rule**: `alignment_metrics`

```python
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
```

### 4. MultiQC Reporting

MultiQC collects all alignment metrics and generates a comprehensive report.

**Rule**: `multiqc_alignment`

```python
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
```

### 5. Alignment Completion Checkpoint

This rule ensures all alignment steps are completed successfully.

**Rule**: `alignment_complete`

```python
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
```

## Input Files

- **Trimmed FASTQ Files**: Generated in the trimming stage
  - `Analysis/Trimmed/{sample}/{sample}_R1_trimmed.fastq.gz`
  - `Analysis/Trimmed/{sample}/{sample}_R2_trimmed.fastq.gz`
- **Trimming Completion Marker**: `Analysis/Trimmed/.trimming_complete`
- **Trimmed Samples Metadata**: `resources/metadata/trimmed_samples.csv`
- **Reference Genome Index**: `resources/genome/star_index`
- **GTF Annotation**: `resources/genome/Homo_sapiens.GRCh38.90.gtf`

### Trimmed Samples Metadata

The `trimmed_samples.csv` file is automatically generated in the trimming stage and contains paths to trimmed FASTQ files. This file is required by the alignment stage to ensure that all trimmed files are properly documented and available.

Example content:
```
sample,R1,R2
sample1,Analysis/Trimmed/sample1/sample1_R1_trimmed.fastq.gz,Analysis/Trimmed/sample1/sample1_R2_trimmed.fastq.gz
sample2,Analysis/Trimmed/sample2/sample2_R1_trimmed.fastq.gz,Analysis/Trimmed/sample2/sample2_R2_trimmed.fastq.gz
```

## Output Files

### STAR Alignment
- `Analysis/Alignment/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam`: Coordinate-sorted BAM file
- `Analysis/Alignment/STAR/{sample}/{sample}.Aligned.toTranscriptome.out.bam`: Transcriptome-aligned BAM file
- `Analysis/Alignment/STAR/{sample}/{sample}.ReadsPerGene.out.tab`: Gene counts
- `Analysis/Alignment/STAR/{sample}/{sample}.SJ.out.tab`: Splice junctions
- `Analysis/Alignment/STAR/{sample}/{sample}.Log.final.out`: STAR summary statistics
- `Analysis/Alignment/STAR/{sample}/{sample}.Log.out`: STAR log file
- `Analysis/Alignment/STAR/{sample}/{sample}.Log.progress.out`: STAR progress log

### BAM Indices
- `Analysis/Alignment/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai`: BAM index

### Alignment Metrics
- `Analysis/Alignment/STAR/{sample}/{sample}.flagstat.txt`: Alignment statistics
- `Analysis/Alignment/STAR/{sample}/{sample}.idxstats.txt`: Index statistics
- `Analysis/Alignment/STAR/{sample}/{sample}.stats.txt`: Detailed alignment statistics

### MultiQC Report
- `Analysis/Alignment/MultiQC/multiqc_report.html`: Aggregated alignment report
- `Analysis/Alignment/MultiQC/multiqc_data/`: Directory containing raw data

### Checkpoint Marker
- `Analysis/Alignment/.alignment_complete`: Marker file indicating alignment completion

## Configuration

The Alignment stage can be configured in `resources/config/params.yaml`:

```yaml
# Reference genome paths
genome_fasta: "resources/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
gtf_file: "resources/genome/Homo_sapiens.GRCh38.90.gtf"

# STAR index parameters
star_index: "resources/genome/star_index"
star_index_threads: 16
star_index_memory: 40000
star_index_time: "04:00:00"
sjdb_overhang: 100  # Read length - 1

# STAR alignment parameters
star_threads: 8
star_memory: 32000
star_align_time: "04:00:00"
star_mismatch: 10
star_multimap: 20
star_score_min: 0.66
star_match_min: 0.66
star_sjdb_overhang_min: 3
star_intron_max: 500000

# BAM processing parameters
bam_index_memory: 4000
bam_index_time: "01:00:00"
metrics_memory: 4000
metrics_time: "01:00:00"

# MultiQC for alignment
multiqc_alignment_memory: 4000
multiqc_alignment_time: "01:00:00"

# Alignment completion
alignment_complete_memory: 1000
alignment_complete_time: "00:10:00"
```

## Running the Alignment Stage

To run only the alignment stage:

```bash
snakemake --cores <N> alignment_stage
```

This will:
1. Run STAR alignment on all samples
2. Index the BAM files
3. Generate alignment metrics
4. Create a MultiQC report
5. Mark alignment as complete

## STAR Index Creation

Before running the alignment, you need to create a STAR index. This can be done using the provided rule:

```bash
snakemake --cores <N> create_star_index
```

This rule requires:
- Reference genome FASTA: `resources/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa`
- GTF annotation: `resources/genome/Homo_sapiens.GRCh38.90.gtf`

## Alignment Parameters

STAR alignment can be customized with various parameters:

- **outFilterMismatchNmax**: Maximum number of mismatches per pair
- **outFilterMultimapNmax**: Maximum number of multiple alignments
- **outFilterScoreMinOverLread**: Minimum alignment score relative to read length
- **outFilterMatchNminOverLread**: Minimum number of matched bases relative to read length
- **alignSJDBoverhangMin**: Minimum overhang for annotated junctions
- **alignIntronMax**: Maximum intron length

## Interpreting Alignment Results

### STAR Log.final.out

The `Log.final.out` file contains summary statistics for each sample, including:
- Number of input reads
- Uniquely mapped reads (%)
- Multi-mapped reads (%)
- Unmapped reads (%)
- Chimeric reads (%)

### Alignment Metrics

The alignment metrics files provide detailed information about the alignments:
- **flagstat**: Summary statistics about the alignments
- **idxstats**: Statistics about the reference sequences
- **stats**: Detailed alignment statistics

### MultiQC Report

The MultiQC report aggregates all alignment metrics into a comprehensive report, making it easy to compare samples and identify outliers.

## Troubleshooting

### Common Issues

1. **STAR Alignment Failures**
   - Check input file paths in the trimmed samples CSV
   - Verify STAR index was created correctly
   - Check memory and thread requirements
   - Review log files for specific errors

2. **BAM Indexing Issues**
   - Verify BAM files were created correctly
   - Check disk space
   - Review log files for errors

3. **Missing Alignment Metrics**
   - Verify BAM files and indices exist
   - Check samtools installation
   - Review log files for errors

### Log Files

- `logs/star/{sample}.log`: STAR alignment logs
- `logs/star/{sample}.index.log`: BAM indexing logs
- `logs/star/{sample}.metrics.log`: Alignment metrics logs
- `logs/star/multiqc.log`: MultiQC log
- `logs/workflow/alignment_complete.log`: Alignment completion log

## Next Steps

After successful completion of the alignment stage, you can proceed with downstream analyses such as:
- Differential expression analysis
- Transcript assembly
- Variant calling
- Pathway analysis

The aligned BAM files and gene counts can be used with tools like:
- DESeq2 or edgeR for differential expression
- StringTie for transcript assembly
- GATK for variant calling
- GSEA for pathway analysis 