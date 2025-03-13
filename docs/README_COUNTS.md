# RNA-seq Feature Counts

This document details the Feature Counts stage of the RNA-seq analysis workflow.

## Overview

The Feature Counts stage quantifies gene expression by counting reads mapped to genomic features. This stage consists of four main steps:

1. **Feature Counts**: Counts reads mapped to genomic features for each sample
2. **Count Merging**: Combines individual count files into a single matrix
3. **Count Normalization**: Generates normalized count values for comparison
4. **MultiQC Reporting**: Aggregates count statistics into a comprehensive report

## Workflow Steps

### 1. Feature Counts

featureCounts quantifies gene expression by counting reads mapped to genomic features defined in a GTF file.

**Rule**: `feature_counts`

```python
rule feature_counts:
    input:
        bam="Analysis/Alignment/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
        bai="Analysis/Alignment/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai",
        alignment_complete="Analysis/Alignment/.main_alignment_complete"
    output:
        counts="Analysis/Counts/FeatureCounts/{sample}/{sample}.counts.txt",
        summary="Analysis/Counts/FeatureCounts/{sample}/{sample}.counts.txt.summary"
    params:
        gtf=config.get("gtf_file", "resources/genome/Homo_sapiens.GRCh38.90.gtf"),
        feature_type=config.get("feature_type", "exon"),
        attribute=config.get("attribute", "gene_id"),
        strand=config.get("strand_specificity", "0"),
        extra_params=config.get("featurecounts_extra", "-p -B -C --fracOverlap 0.2")
    threads: config.get("featurecounts_threads", 4)
    resources:
        mem_mb=config.get("featurecounts_memory", 8000)
```

This rule:
- Takes aligned BAM files as input
- Counts reads mapped to features (default: exons)
- Groups counts by attribute (default: gene_id)
- Handles paired-end reads with proper fragment counting
- Produces a count file and summary statistics for each sample

### 2. Count Merging

This step combines individual count files into a single matrix for easier analysis.

**Rule**: `merge_counts`

```python
rule merge_counts:
    input:
        counts=expand("Analysis/Counts/FeatureCounts/{sample}/{sample}.counts.txt", 
                     sample=SAMPLES.keys())
    output:
        merged="Analysis/Counts/FeatureCounts/merged_gene_counts.txt",
        normalized="Analysis/Counts/FeatureCounts/merged_gene_counts.normalized.txt"
    resources:
        mem_mb=config.get("merge_counts_memory", 8000)
    script:
        "scripts/merge_counts.py"
```

This rule:
- Takes all individual count files as input
- Creates a merged count matrix with genes as rows and samples as columns
- Generates a normalized count matrix (CPM normalization)
- Uses a Python script for flexible merging and normalization

### 3. MultiQC Reporting

MultiQC collects all count summaries and generates a comprehensive report.

**Rule**: `multiqc_counts`

```python
rule multiqc_counts:
    input:
        summaries=expand("Analysis/Counts/FeatureCounts/{sample}/{sample}.counts.txt.summary", 
                        sample=SAMPLES.keys())
    output:
        html="Analysis/Counts/MultiQC/multiqc_report.html",
        data_dir=directory("Analysis/Counts/MultiQC/multiqc_data")
    params:
        outdir="Analysis/Counts/MultiQC"
    resources:
        mem_mb=config.get("multiqc_counts_memory", 4000)
```

This rule:
- Takes all count summary files as input
- Generates a comprehensive report of counting statistics
- Provides visualizations for comparing samples

### 4. Counts Completion Checkpoint

This rule ensures all counting steps are completed successfully.

**Rule**: `counts_complete`

```python
rule counts_complete:
    input:
        counts=expand("Analysis/Counts/FeatureCounts/{sample}/{sample}.counts.txt", 
                     sample=SAMPLES.keys()),
        merged="Analysis/Counts/FeatureCounts/merged_gene_counts.txt",
        normalized="Analysis/Counts/FeatureCounts/merged_gene_counts.normalized.txt",
        multiqc="Analysis/Counts/MultiQC/multiqc_report.html"
    output:
        touch("Analysis/Counts/.counts_complete")
```

This rule:
- Verifies all count files were generated
- Confirms the merged and normalized matrices exist
- Creates a marker file indicating completion

## Input Files

- **Aligned BAM Files**: Generated in the alignment stage
  - `Analysis/Alignment/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam`
- **BAM Indices**: Generated in the alignment stage
  - `Analysis/Alignment/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai`
- **Alignment Completion Marker**: `Analysis/Alignment/.main_alignment_complete`
- **GTF Annotation**: `resources/genome/Homo_sapiens.GRCh38.90.gtf`

## Output Files

### Feature Counts
- `Analysis/Counts/FeatureCounts/{sample}/{sample}.counts.txt`: Count file for each sample
- `Analysis/Counts/FeatureCounts/{sample}/{sample}.counts.txt.summary`: Summary statistics

### Merged Counts
- `Analysis/Counts/FeatureCounts/merged_gene_counts.txt`: Merged count matrix
- `Analysis/Counts/FeatureCounts/merged_gene_counts.normalized.txt`: Normalized count matrix

### MultiQC Report
- `Analysis/Counts/MultiQC/multiqc_report.html`: Aggregated count report
- `Analysis/Counts/MultiQC/multiqc_data/`: Directory containing raw data

### Checkpoint Marker
- `Analysis/Counts/.counts_complete`: Marker file indicating counts completion

## Configuration

The Feature Counts stage can be configured in `resources/config/params.yaml`:

```yaml
# Feature counts parameters
featurecounts_threads: 4
featurecounts_memory: 8000
merge_counts_memory: 4000
multiqc_counts_memory: 4000

# Feature counts options
feature_type: "exon"          # Feature type to count (exon, gene, transcript)
attribute: "gene_id"          # Attribute to use for grouping features
strand_specificity: "0"       # Strand specificity: 0 (unstranded), 1 (stranded), 2 (reversely stranded)
featurecounts_extra: "-p -B -C --fracOverlap 0.2"  # Additional parameters for featureCounts
```

## Detailed Parameter Explanations

### featureCounts Resource Parameters

featureCounts is a highly efficient general-purpose read summarization program that counts mapped reads for genomic features.

| Parameter | Description | Default |
|-----------|-------------|---------|
| `featurecounts_threads` | Number of CPU threads to use for featureCounts. Higher values can significantly speed up processing, especially for large BAM files. | 4 |
| `featurecounts_memory` | Memory allocation in MB for featureCounts. Large BAM files or complex annotations may require more memory. | 8000 |
| `merge_counts_memory` | Memory allocation in MB for merging count files. This increases with the number of samples and genes. | 4000 |
| `multiqc_counts_memory` | Memory allocation in MB for MultiQC when processing count reports. | 4000 |

### featureCounts Feature Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `feature_type` | Type of feature to count. This should match the feature type in your GTF file. Common values are "exon", "gene", or "transcript". | "exon" |
| `attribute` | Attribute in the GTF file to use for grouping features. This determines how reads are assigned to genes. Common values are "gene_id", "gene_name", or "transcript_id". | "gene_id" |
| `strand_specificity` | Strand-specific counting mode:<br>- 0: Unstranded (counts reads on both strands)<br>- 1: Stranded (counts reads on the same strand)<br>- 2: Reversely stranded (counts reads on the opposite strand)<br>This should match your library preparation protocol. | "0" |

### featureCounts Advanced Parameters

The `featurecounts_extra` parameter allows you to specify additional options for featureCounts. The default value is:

```
-p -B -C --fracOverlap 0.2
```

These options mean:

| Option | Description |
|--------|-------------|
| `-p` | Count fragments (paired-end) instead of individual reads. This is recommended for paired-end data to avoid counting both reads of a pair. |
| `-B` | Only count read pairs that have both ends mapped. This ensures higher quality counts by excluding partially mapped fragments. |
| `-C` | Don't count reads that map to multiple features. This reduces ambiguity in the counts. |
| `--fracOverlap 0.2` | Minimum fraction of overlap required between a read and a feature. A read must overlap at least 20% of a feature to be counted. |

### Other Useful featureCounts Options

You can add these to the `featurecounts_extra` parameter as needed:

| Option | Description |
|--------|-------------|
| `-O` | Allow reads to be assigned to multiple features (not recommended for standard RNA-seq). |
| `-M` | Count multi-mapping reads. By default, only uniquely mapped reads are counted. |
| `--primary` | Count primary alignments only (exclude secondary and supplementary alignments). |
| `--ignoreDup` | Ignore duplicate reads (those marked with the 0x400 flag). |
| `--minOverlap <N>` | Minimum number of overlapping bases required between a read and a feature. |
| `--countReadPairs` | Count read pairs instead of reads (similar to -p but with different behavior). |
| `--donotsort` | Do not sort BAM files by read name (saves time if BAMs are already name-sorted). |

### Count Merging Parameters

The count merging step combines individual count files into a single matrix and generates normalized counts.

| Parameter | Description | Default |
|-----------|-------------|---------|
| `merge_counts_memory` | Memory allocation in MB for merging count files. | 4000 |

The merging script performs the following operations:
1. Extracts gene counts from individual featureCounts output files
2. Combines them into a single matrix with genes as rows and samples as columns
3. Generates a normalized count matrix using CPM (Counts Per Million) normalization
4. Handles metadata and sample information from the input files

### MultiQC Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `multiqc_counts_memory` | Memory allocation in MB for MultiQC when processing count reports. | 4000 |

## Running the Feature Counts Stage

To run only the feature counts stage:

```bash
snakemake --cores <N> count_all
```

This will:
1. Run featureCounts on all samples
2. Merge individual count files
3. Generate normalized counts
4. Create a MultiQC report
5. Mark counts as complete

## Feature Counts Parameters

featureCounts can be customized with various parameters:

- **feature_type**: Type of feature to count (exon, gene, transcript)
- **attribute**: Attribute to use for grouping features (gene_id, gene_name)
- **strand_specificity**: Strand-specific counting mode
  - 0: Unstranded (default)
  - 1: Stranded
  - 2: Reversely stranded
- **extra_params**: Additional parameters for featureCounts
  - `-p`: Count fragments (paired-end) instead of reads
  - `-B`: Only count read pairs that have both ends mapped
  - `-C`: Don't count reads that map to multiple features
  - `--fracOverlap 0.2`: Minimum fraction of overlap required

## Interpreting Count Results

### Count Files

Each sample's count file contains:
- Genomic feature information (chromosome, start, end, strand)
- Feature annotation (gene ID, gene name)
- Raw count values

### Summary Statistics

The summary files provide information about:
- Total reads processed
- Successfully assigned reads
- Unassigned reads (with reasons)
  - Ambiguous
  - No features
  - Low mapping quality
  - Chimeric reads
  - Multimapping

### Merged Count Matrix

The merged count matrix contains:
- Rows: Genes or features
- Columns: Samples
- Values: Raw count values

### Normalized Count Matrix

The normalized count matrix contains:
- Rows: Genes or features
- Columns: Samples
- Values: Counts per million (CPM) normalized values

## Troubleshooting

### Common Issues

1. **Feature Counts Failures**
   - Check input BAM files
   - Verify GTF file format and compatibility
   - Check memory requirements
   - Review log files for specific errors

2. **Low Assignment Rate**
   - Check strand specificity setting
   - Verify feature type matches experimental design
   - Check GTF annotation version matches reference genome

3. **Merging Issues**
   - Verify all count files were generated
   - Check disk space
   - Review log files for errors

### Log Files

- `logs/counts/{sample}.featurecounts.log`: Feature counts logs for each sample
- `logs/counts/merge_counts.log`: Count merging log
- `logs/counts/multiqc.log`: MultiQC log
- `logs/workflow/counts_complete.log`: Counts completion log

## Next Steps

After successful completion of the feature counts stage, you can proceed with downstream analyses such as:
- Differential expression analysis with DESeq2 or edgeR
- Pathway analysis
- Gene set enrichment analysis
- Visualization and exploration of expression patterns

The merged count matrix is ready for import into R or Python for these analyses. 