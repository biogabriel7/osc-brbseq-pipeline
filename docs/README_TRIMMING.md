# RNA-seq Trimming

This document details the Trimming stage of the RNA-seq analysis workflow.

## Overview

The Trimming stage processes raw sequencing reads to remove adapters, low-quality bases, and other artifacts that could interfere with downstream analysis. This stage consists of three main steps:

1. **fastp Trimming**: Performs adapter and quality trimming with sample-specific parameters
2. **MultiQC Reporting**: Aggregates trimming reports into a comprehensive summary
3. **Trimmed Samples Metadata Generation**: Automatically creates a CSV file with paths to trimmed files

## Workflow Steps

### 1. fastp Trimming

fastp performs adapter removal and quality trimming using parameters determined in the QC stage.

**Rule**: `fastp_with_dependency`

```python
rule fastp_with_dependency:
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
        fastp_opts=get_fastp_params,
        time=config.get("trimming_time", "02:00:00")
    threads: config.get("trimming_threads", 4)
    resources:
        mem_mb=config.get("trimming_memory", 8000)
```

The `get_fastp_params` function retrieves sample-specific parameters from the `trimming_params.json` file generated in the QC stage.

### 2. MultiQC Reporting

MultiQC collects all fastp reports and generates a comprehensive summary.

**Rule**: `multiqc_fastp`

```python
rule multiqc_fastp:
    input:
        fastp_reports=expand("Analysis/QC/Trimming/Reports/{sample}_fastp.json",
                           sample=SAMPLES.keys())
    output:
        html="Analysis/QC/Trimming/MultiQC/multiqc_report.html",
        data_dir=directory("Analysis/QC/Trimming/MultiQC/multiqc_data")
    log:
        "logs/trimming/multiqc_fastp.log"
    params:
        time=config.get("multiqc_fastp_time", "00:30:00")
    resources:
        mem_mb=config.get("multiqc_fastp_memory", 4000)
```

### 3. Trimming Completion Checkpoint

This rule ensures all trimming steps are completed successfully.

**Rule**: `trimming_complete`

```python
rule trimming_complete:
    input:
        trimmed=expand("Analysis/Trimmed/{sample}/{sample}_R1_trimmed.fastq.gz",
                      sample=SAMPLES.keys()),
        multiqc="Analysis/QC/Trimming/MultiQC/multiqc_report.html"
    output:
        touch("Analysis/Trimmed/.trimming_complete")
    log:
        "logs/workflow/trimming_complete.log"
    params:
        time=config.get("trimming_complete_time", "00:10:00")
    resources:
        mem_mb=config.get("trimming_complete_memory", 1000)
```

### 4. Trimmed Samples Metadata Generation

This rule automatically generates a CSV file with paths to trimmed FASTQ files, which is used by the alignment stage.

**Rule**: `generate_trimmed_samples_csv`

```python
rule generate_trimmed_samples_csv:
    input:
        trimming_complete="Analysis/Trimmed/.trimming_complete"
    output:
        csv="resources/metadata/trimmed_samples.csv"
    log:
        "logs/workflow/generate_trimmed_samples.log"
    params:
        time=config.get("generate_trimmed_csv_time", "00:10:00")
    resources:
        mem_mb=config.get("generate_trimmed_csv_memory", 1000)
```

This rule:
1. Finds all trimmed R1 files in the `Analysis/Trimmed` directory
2. Extracts the sample name from the path
3. Finds the corresponding R2 file
4. Creates a CSV file with sample name, R1 path, and R2 path

## Input Files

- **Raw FASTQ Files**: Specified in the metasheet (`resources/config/metasheet.csv`)
- **QC Completion Marker**: `Analysis/QC/.qc_complete`
- **Trimming Parameters**: `Analysis/QC/Trimming/trimming_params.json`

## Output Files

### Trimmed FASTQ Files
- `Analysis/Trimmed/{sample}/{sample}_R1_trimmed.fastq.gz`: Trimmed R1 reads
- `Analysis/Trimmed/{sample}/{sample}_R2_trimmed.fastq.gz`: Trimmed R2 reads

### Trimming Reports
- `Analysis/QC/Trimming/Reports/{sample}_fastp.html`: HTML report for each sample
- `Analysis/QC/Trimming/Reports/{sample}_fastp.json`: JSON report for each sample

### MultiQC Report
- `Analysis/QC/Trimming/MultiQC/multiqc_report.html`: Aggregated trimming report
- `Analysis/QC/Trimming/MultiQC/multiqc_data/`: Directory containing raw data

### Checkpoint Markers
- `Analysis/Trimmed/.trimming_complete`: Marker file indicating trimming completion

### Trimmed Samples Metadata
- `resources/metadata/trimmed_samples.csv`: CSV file with paths to trimmed FASTQ files

Example content of `trimmed_samples.csv`:
```
sample,R1,R2
sample1,Analysis/Trimmed/sample1/sample1_R1_trimmed.fastq.gz,Analysis/Trimmed/sample1/sample1_R2_trimmed.fastq.gz
sample2,Analysis/Trimmed/sample2/sample2_R1_trimmed.fastq.gz,Analysis/Trimmed/sample2/sample2_R2_trimmed.fastq.gz
```

## Configuration

The Trimming stage can be configured in `resources/config/params.yaml`:

```yaml
# Trimming parameters
trimming_threads: 4
trimming_memory: 8000
trimming_time: "02:00:00"

# Default trimming parameters (used if not specified in trimming_params.json)
leading_quality: 3
trailing_quality: 3
sliding_window: "4:15"
min_length: 36

# MultiQC for trimming
multiqc_fastp_memory: 4000
multiqc_fastp_time: "00:30:00"

# Trimming completion
trimming_complete_memory: 1000
trimming_complete_time: "00:10:00"

# Trimmed samples CSV generation
generate_trimmed_csv_memory: 1000
generate_trimmed_csv_time: "00:10:00"
```

## Running the Trimming Stage

To run only the trimming stage:

```bash
snakemake --cores <N> trimming_stage
```

This will:
1. Run fastp on all samples
2. Generate MultiQC report
3. Create the trimming completion marker
4. Generate the trimmed samples CSV file

## Trimming Parameters

The `get_fastp_params` function retrieves sample-specific parameters from the `trimming_params.json` file. If parameters for a sample are not found, default values from `params.yaml` are used.

Parameters include:
- **leading_quality**: Quality threshold for trimming bases from the start of reads
- **trailing_quality**: Quality threshold for trimming bases from the end of reads
- **min_length**: Minimum read length to keep after trimming
- **sliding_window**: Window size and quality threshold for sliding window trimming
- **deduplication**: Whether to remove duplicate reads
- **trim_poly_g**: Whether to trim poly-G tails (common in Illumina NextSeq/NovaSeq data)

## Automatic Trimmed Samples Metadata Generation

The `generate_trimmed_samples_csv` rule automatically creates a CSV file with paths to trimmed FASTQ files. This file is used by the alignment stage to locate the trimmed reads.

The rule:
1. Runs after all trimming is complete
2. Finds all trimmed R1 files in the `Analysis/Trimmed` directory
3. Extracts the sample name from the path
4. Finds the corresponding R2 file
5. Creates a CSV file with sample name, R1 path, and R2 path

This automation ensures that:
- All trimmed files are properly documented
- The alignment stage has accurate paths to trimmed files
- No manual intervention is required between stages

## Troubleshooting

### Common Issues

1. **Trimming Failures**
   - Check input file paths in metasheet
   - Verify fastp is installed and in PATH
   - Check log files for specific errors

2. **Missing Trimmed Files**
   - Verify fastp completed successfully
   - Check disk space
   - Review log files for errors

3. **Trimmed Samples CSV Generation Issues**
   - Verify all trimming jobs completed
   - Check directory structure
   - Review log file for specific errors

### Log Files

- `logs/trimming/{sample}.log`: fastp logs for each sample
- `logs/trimming/multiqc_fastp.log`: MultiQC log
- `logs/workflow/trimming_complete.log`: Trimming completion log
- `logs/workflow/generate_trimmed_samples.log`: Trimmed samples CSV generation log

## Next Steps

After successful completion of the trimming stage, the workflow will proceed to the alignment stage, which uses the trimmed reads and the automatically generated `trimmed_samples.csv` file to align reads to the reference genome. 