# RNA-seq Quality Control

This document details the Quality Control (QC) stage of the RNA-seq analysis workflow.

## Overview

The QC stage performs quality assessment of raw sequencing data and automatically determines optimal trimming parameters for each sample. This stage consists of three main steps:

1. **FastQC Analysis**: Evaluates the quality of raw sequencing reads
2. **MultiQC Reporting**: Aggregates FastQC results into a comprehensive report
3. **Parameter Generation**: Analyzes QC metrics to determine sample-specific trimming parameters

## Workflow Steps

### 1. FastQC Analysis

FastQC examines raw FASTQ files to identify quality issues such as:
- Base quality scores
- Per-sequence quality
- GC content distribution
- Sequence duplication levels
- Adapter content
- Overrepresented sequences

**Rule**: `fastqc`

```python
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
        outdir="Analysis/QC/FastQC/{sample}",
        time=config.get("fastqc_time", "01:00:00")
    threads: config.get("fastqc_threads", 2)
    resources:
        mem_mb=config.get("fastqc_memory", 4000)
```

### 2. MultiQC Reporting

MultiQC collects all FastQC results and generates a comprehensive report for easy comparison across samples.

**Rule**: `multiqc`

```python
rule multiqc:
    input:
        fastqc_outputs=expand("Analysis/QC/FastQC/{sample}/{sample}_{read}_fastqc.zip",
                             sample=SAMPLES.keys(),
                             read=["R1", "R2"])
    output:
        html="Analysis/QC/MultiQC/multiqc_report.html",
        data_dir=directory("Analysis/QC/MultiQC/multiqc_data")
    log:
        "logs/multiqc/multiqc.log"
    params:
        outdir="Analysis/QC/MultiQC",
        time=config.get("multiqc_time", "00:30:00")
    resources:
        mem_mb=config.get("multiqc_memory", 4000)
```

### 3. Parameter Generation

This step analyzes the FastQC results to determine optimal trimming parameters for each sample, which are stored in a JSON file.

**Rule**: `generate_trimming_params`

```python
rule generate_trimming_params:
    input:
        multiqc_html="Analysis/QC/MultiQC/multiqc_report.html",
        config=config.get("params_path", "resources/config/params.yaml")
    output:
        params="Analysis/QC/Trimming/trimming_params.json"
    log:
        "logs/trimming/generate_params.log"
    params:
        time=config.get("params_time", "00:30:00")
```

### 4. QC Completion Checkpoint

This rule ensures all QC steps are completed successfully before proceeding to the trimming stage.

**Rule**: `qc_complete`

```python
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
```

## Input Files

- **Raw FASTQ Files**: Specified in the metasheet (`resources/config/metasheet.csv`)
  ```
  sample,R1,R2
  sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz
  ```

## Output Files

### FastQC Reports
- `Analysis/QC/FastQC/{sample}/{sample}_R1_fastqc.html`: HTML report for R1 reads
- `Analysis/QC/FastQC/{sample}/{sample}_R1_fastqc.zip`: Raw data for R1 reads
- `Analysis/QC/FastQC/{sample}/{sample}_R2_fastqc.html`: HTML report for R2 reads
- `Analysis/QC/FastQC/{sample}/{sample}_R2_fastqc.zip`: Raw data for R2 reads

### MultiQC Report
- `Analysis/QC/MultiQC/multiqc_report.html`: Aggregated QC report
- `Analysis/QC/MultiQC/multiqc_data/`: Directory containing raw data

### Trimming Parameters
- `Analysis/QC/Trimming/trimming_params.json`: Sample-specific trimming parameters

### Checkpoint Marker
- `Analysis/QC/.qc_complete`: Marker file indicating QC completion

## Configuration

The QC stage can be configured in `resources/config/params.yaml`:

```yaml
# FastQC specific parameters
fastqc_threads: 2
fastqc_memory: 4000
fastqc_time: "01:00:00"

# MultiQC parameters
multiqc_memory: 4000
multiqc_time: "00:30:00"

# Parameter generation
params_time: "00:30:00"

# QC completion
qc_complete_memory: 1000
qc_complete_time: "00:10:00"
```

## Detailed Parameter Explanations

### FastQC Parameters

FastQC is a quality control tool for high throughput sequence data. It runs a set of analyses on your FASTQ files and provides a comprehensive report.

| Parameter | Description | Default |
|-----------|-------------|---------|
| `fastqc_threads` | Number of CPU threads to use for FastQC analysis. Increasing this value can speed up processing for large files. | 2 |
| `fastqc_memory` | Memory allocation in MB for FastQC. Most samples require minimal memory, but larger files may need more. | 4000 |
| `fastqc_time` | Maximum runtime for FastQC jobs in HH:MM:SS format. | "01:00:00" |

FastQC performs several quality checks on raw sequence data:
- **Basic Statistics**: Summary statistics about the file
- **Per Base Sequence Quality**: Quality scores across all bases at each position
- **Per Sequence Quality Scores**: Quality score distribution over all sequences
- **Per Base Sequence Content**: Proportion of each base at each position
- **Per Sequence GC Content**: GC content across all sequences
- **Per Base N Content**: Proportion of N bases at each position
- **Sequence Length Distribution**: Distribution of sequence lengths
- **Sequence Duplication Levels**: Degree of duplication for sequences
- **Overrepresented Sequences**: Sequences that appear more than expected
- **Adapter Content**: Adapter sequences present in the data
- **Kmer Content**: Sequences that are overrepresented at specific positions

### MultiQC Parameters

MultiQC aggregates results from multiple bioinformatics analyses into a single report.

| Parameter | Description | Default |
|-----------|-------------|---------|
| `multiqc_memory` | Memory allocation in MB for MultiQC. | 4000 |
| `multiqc_time` | Maximum runtime for MultiQC jobs in HH:MM:SS format. | "00:30:00" |

### Parameter Generation Parameters

The parameter generation step analyzes FastQC results to determine optimal trimming parameters for each sample.

| Parameter | Description | Default |
|-----------|-------------|---------|
| `params_time` | Maximum runtime for parameter generation in HH:MM:SS format. | "00:30:00" |

The parameter generation script analyzes the following metrics to determine trimming parameters:
- **Adapter Content**: Determines if adapter trimming is needed
- **Per Base Sequence Quality**: Determines quality thresholds for trimming
- **Overrepresented Sequences**: Identifies potential contaminants
- **GC Content**: Identifies potential bias or contamination
- **Sequence Duplication**: Determines if deduplication is needed
- **Poly-G Content**: Determines if poly-G trimming is needed (common in NextSeq/NovaSeq data)

The script generates a JSON file with sample-specific parameters that will be used in the trimming stage:
- `leading_quality`: Quality threshold for trimming bases from the start of reads
- `trailing_quality`: Quality threshold for trimming bases from the end of reads
- `min_length`: Minimum read length to keep after trimming
- `sliding_window`: Window size and quality threshold for sliding window trimming
- `deduplication`: Whether to remove duplicate reads
- `trim_poly_g`: Whether to trim poly-G tails

## Running the QC Stage

To run only the QC stage:

```bash
snakemake --cores <N> qc_stage
```

## Interpreting QC Results

### FastQC Reports

FastQC reports provide detailed quality metrics for each sample. Key sections to review:

1. **Per Base Sequence Quality**: Look for high quality scores (green area)
2. **Per Sequence Quality Scores**: Should show a peak at high quality
3. **Per Base Sequence Content**: Should show balanced nucleotide distribution
4. **Adapter Content**: Identifies adapter contamination
5. **Overrepresented Sequences**: Identifies potential contamination or bias

### MultiQC Summary

The MultiQC report provides a comparative view of all samples, making it easier to identify:
- Outlier samples with quality issues
- Systematic problems affecting multiple samples
- Patterns in quality metrics across the dataset

### Trimming Parameters

The `trimming_params.json` file contains sample-specific parameters determined from QC results:

```json
{
  "sample1": {
    "leading_quality": 3,
    "trailing_quality": 3,
    "min_length": 36,
    "sliding_window": "4:15",
    "deduplication": false,
    "trim_poly_g": true
  },
  "sample2": {
    ...
  }
}
```

These parameters will be used in the trimming stage to optimize read quality for each sample.

## Troubleshooting

### Common Issues

1. **Missing FastQC Reports**
   - Check input file paths in metasheet
   - Verify FastQC is installed and in PATH
   - Check log files for errors

2. **Poor Quality Metrics**
   - High adapter content: Will be addressed in trimming
   - Low quality scores: May require more aggressive trimming
   - Overrepresented sequences: Check for contamination

3. **Parameter Generation Failures**
   - Verify MultiQC completed successfully
   - Check MultiQC data directory for required metrics files
   - Review log file for specific errors

### Log Files

- `logs/fastqc/{sample}.log`: FastQC logs for each sample
- `logs/multiqc/multiqc.log`: MultiQC log
- `logs/trimming/generate_params.log`: Parameter generation log
- `logs/workflow/qc_complete.log`: QC completion log

## Next Steps

After successful completion of the QC stage, the workflow will proceed to the trimming stage, which uses the generated parameters to optimize read quality. 