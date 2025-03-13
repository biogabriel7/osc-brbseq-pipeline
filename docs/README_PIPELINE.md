# RNA-seq Analysis Pipeline

## Overview

This RNA-seq analysis pipeline provides a comprehensive workflow for processing and analyzing RNA sequencing data. The pipeline is implemented using Snakemake, a workflow management system that ensures reproducibility and scalability.

## Pipeline Stages

The pipeline consists of five main stages:

1. **[Quality Control (QC)](README_QC.md)**: Assesses the quality of raw sequencing data and determines optimal trimming parameters
2. **[Trimming](README_TRIMMING.md)**: Removes adapters and low-quality bases from raw reads
3. **[Alignment](README_ALIGNMENT.md)**: Maps trimmed reads to a reference genome
4. **[Feature Counts](README_COUNTS.md)**: Quantifies gene expression by counting reads mapped to genomic features
5. **Workflow Completion**: Ensures all stages are completed successfully

Each stage has its own checkpoint to ensure that all required outputs are generated before proceeding to the next stage.

## Directory Structure

The pipeline creates the following directory structure:

```
project_root/
├── Analysis/
│   ├── QC/
│   │   ├── FastQC/
│   │   ├── MultiQC/
│   │   └── Trimming/
│   ├── Trimmed/
│   ├── Alignment/
│   │   ├── STAR/
│   │   └── MultiQC/
│   └── Counts/
│       ├── FeatureCounts/
│       └── MultiQC/
├── logs/
│   ├── fastqc/
│   ├── trimming/
│   ├── star/
│   ├── counts/
│   └── workflow/
├── resources/
│   ├── config/
│   │   ├── params.yaml
│   │   ├── metasheet.csv
│   │   └── cluster.yaml
│   ├── genome/
│   └── metadata/
├── workflow/
│   ├── qc_params.snakefile
│   ├── fastp_trimming.snakefile
│   ├── star_alignment.snakefile
│   ├── feature_counts.snakefile
│   ├── scripts/
│   └── common/
└── Snakefile
```

## Installation and Setup

### Requirements

- Snakemake (≥ 6.0)
- FastQC
- MultiQC
- fastp
- STAR
- Samtools
- Subread (for featureCounts)
- Python (≥ 3.6)

### Setup Instructions

1. Clone the repository:
   ```bash
   git clone https://github.com/username/rnaseq-pipeline.git
   cd rnaseq-pipeline
   ```

2. Create a conda environment:
   ```bash
   conda env create -f environment.yaml
   conda activate rnaseq
   ```

3. Configure the pipeline:
   - Edit `resources/config/params.yaml` to set pipeline parameters
   - Create or update `resources/config/metasheet.csv` with sample information

4. Prepare reference files:
   - Download reference genome and annotation files
   - Create STAR index (or use the provided rule)

## Running the Pipeline

### Basic Usage

To run the entire pipeline:

```bash
snakemake --cores <N> --use-conda
```

### Running Specific Stages

To run specific stages of the pipeline:

```bash
# Run only QC
snakemake --cores <N> qc_all

# Run through trimming
snakemake --cores <N> trim_all

# Run through alignment
snakemake --cores <N> align_all

# Run through feature counts
snakemake --cores <N> count_all
```

### Cluster Execution

To run the pipeline on a cluster (e.g., SLURM):

```bash
snakemake --profile slurm --cluster-config resources/config/cluster.yaml
```

## Configuration

### Sample Configuration

Samples are defined in `resources/config/metasheet.csv`:

```csv
sample,R1,R2
sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz
sample2,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz
```

### Pipeline Parameters

Pipeline parameters are defined in `resources/config/params.yaml`:

```yaml
# Resource parameters
fastqc_threads: 2
fastqc_memory: 4000
trimming_threads: 4
trimming_memory: 32000
star_threads: 8
star_memory: 32000
featurecounts_threads: 4
featurecounts_memory: 8000

# Reference files
genome_fasta: "resources/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
gtf_file: "resources/genome/Homo_sapiens.GRCh38.90.gtf"
star_index: "resources/genome/star_index"

# STAR parameters
sjdb_overhang: 100
star_mismatch: 10
star_multimap: 20

# Feature counts parameters
feature_type: "exon"
attribute: "gene_id"
strand_specificity: "0"
```

### Cluster Configuration

Cluster parameters are defined in `resources/config/cluster.yaml`:

```yaml
__default__:
  partition: general
  time: "04:00:00"
  nodes: 1
  ntasks: 1
  mem: "4G"

star_align:
  partition: bigmem
  time: "08:00:00"
  mem: "32G"
  ntasks: 8
```

## Pipeline Outputs

### Quality Control
- FastQC reports for each sample
- MultiQC summary report
- Trimming parameters for each sample

### Trimming
- Trimmed FASTQ files
- Trimming reports
- MultiQC summary report

### Alignment
- Aligned BAM files
- BAM indices
- Alignment statistics
- MultiQC summary report

### Feature Counts
- Count files for each sample
- Merged count matrix
- Normalized count matrix
- MultiQC summary report

## Workflow Checkpoints

The pipeline includes checkpoint markers to ensure each stage completes successfully:

- `Analysis/QC/.qc_complete`: QC stage completed
- `Analysis/Trimmed/.trimming_complete`: Trimming stage completed
- `Analysis/Alignment/.main_alignment_complete`: Alignment stage completed
- `Analysis/Counts/.counts_complete`: Feature counts stage completed
- `Analysis/.workflow_complete`: Entire workflow completed

## Troubleshooting

### Common Issues

1. **Missing Input Files**
   - Check paths in metasheet.csv
   - Verify file permissions

2. **Resource Limitations**
   - Increase memory or threads in params.yaml
   - Use cluster execution for large datasets

3. **Software Compatibility**
   - Ensure all tools are installed and in PATH
   - Check version compatibility

### Log Files

Log files are stored in the `logs/` directory, organized by pipeline stage:

- `logs/fastqc/`: FastQC logs
- `logs/trimming/`: Trimming logs
- `logs/star/`: Alignment logs
- `logs/counts/`: Feature counts logs
- `logs/workflow/`: Workflow checkpoint logs

## Advanced Usage

### Test Mode

The pipeline supports a test mode for quick validation:

```bash
snakemake --config test_mode=True --cores <N>
```

### Custom Parameters

You can override parameters at runtime:

```bash
snakemake --config star_threads=16 star_memory=64000 --cores <N>
```

### Resuming Failed Runs

To resume a failed run:

```bash
snakemake --cores <N> --rerun-incomplete
```

## Stage-Specific Documentation

For detailed information about each stage, see the following documents:

- [Quality Control Documentation](README_QC.md)
- [Trimming Documentation](README_TRIMMING.md)
- [Alignment Documentation](README_ALIGNMENT.md)
- [Feature Counts Documentation](README_COUNTS.md) 