# RNA-seq Pipeline Workflow

This directory contains the Snakemake workflow files that define the RNA-seq analysis pipeline.

## Workflow Structure

The workflow is organized into modular Snakemake files, each handling a specific stage of the analysis:

1. **qc_params.snakefile**: Quality control and parameter generation
2. **fastp_trimming.snakefile**: Read trimming and filtering
3. **star_alignment.snakefile**: Genome alignment
4. **feature_counts.snakefile**: Gene expression quantification

## Main Workflow Files

### qc_params.snakefile

This file contains rules for:
- Running FastQC on raw reads
- Generating MultiQC reports
- Analyzing QC metrics to determine optimal trimming parameters

Key rules:
- `fastqc`: Runs FastQC on raw FASTQ files
- `multiqc`: Aggregates FastQC results
- `generate_trimming_params`: Determines sample-specific trimming parameters

### fastp_trimming.snakefile

This file contains rules for:
- Trimming adapters and low-quality bases
- Generating trimming reports
- Creating MultiQC reports for trimming results

Key rules:
- `fastp_trim`: Performs adapter and quality trimming
- `multiqc_fastp`: Aggregates trimming reports
- `generate_trimmed_samples_csv`: Creates metadata for trimmed files

### star_alignment.snakefile

This file contains rules for:
- Aligning trimmed reads to the reference genome
- Indexing BAM files
- Generating alignment metrics
- Creating MultiQC reports for alignment results

Key rules:
- `create_star_index`: Creates STAR genome index
- `star_align`: Aligns reads to the genome
- `index_bam`: Indexes BAM files
- `alignment_metrics`: Generates alignment statistics
- `multiqc_alignment`: Aggregates alignment reports

### feature_counts.snakefile

This file contains rules for:
- Counting reads mapped to genomic features
- Merging count files into a matrix
- Normalizing counts
- Creating MultiQC reports for count results

Key rules:
- `feature_counts`: Counts reads mapped to features
- `merge_counts`: Combines individual count files
- `multiqc_counts`: Aggregates count reports
- `counts_complete`: Marks completion of counting stage

## Utility Directories

- **scripts/**: Contains Python scripts used by the workflow
- **common/**: Contains shared utility functions and modules
- **utils/**: Contains additional utility scripts

## Workflow Execution

The main Snakefile in the root directory imports these modular files and defines the overall workflow execution order. Each stage has a checkpoint rule that ensures all required outputs are generated before proceeding to the next stage.

## Extending the Workflow

To add new functionality:
1. Create a new Snakefile in the workflow directory
2. Define rules for the new functionality
3. Include the new Snakefile in the main Snakefile
4. Update the workflow order in the main Snakefile

## Dependencies Between Stages

The workflow enforces dependencies between stages using checkpoint files:
- QC must complete before trimming starts
- Trimming must complete before alignment starts
- Alignment must complete before feature counts starts

This ensures that each stage has all the required inputs from the previous stage. 