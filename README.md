# RNA-seq Analysis Workflow

A comprehensive Snakemake workflow for RNA-seq data analysis, from quality control to alignment.

## Overview

This workflow automates the RNA-seq analysis process with the following stages:

1. **Quality Control**: FastQC analysis of raw reads and MultiQC reporting
2. **Parameter Generation**: Automatic determination of optimal trimming parameters based on QC results
3. **Trimming**: Adapter and quality trimming using fastp with sample-specific parameters
4. **Alignment**: STAR alignment of trimmed reads to reference genome
5. **Post-alignment Analysis**: BAM indexing, metrics collection, and MultiQC reporting

The workflow is designed to be modular, with clear dependencies between stages and automatic generation of necessary metadata files.

## Directory Structure

```
.
├── Analysis/                  # Main output directory
│   ├── QC/                    # Quality control outputs
│   │   ├── FastQC/           # FastQC reports
│   │   ├── MultiQC/          # MultiQC reports for raw data
│   │   └── Trimming/         # Trimming parameters and reports
│   ├── Trimmed/              # Trimmed FASTQ files
│   └── Alignment/            # Alignment outputs
│       ├── STAR/             # STAR alignment files
│       └── MultiQC/          # MultiQC reports for alignment
├── logs/                      # Log files for all steps
├── resources/                 # Input resources
│   ├── config/               # Configuration files
│   │   ├── params.yaml       # Pipeline parameters
│   │   └── metasheet.csv     # Sample metadata
│   ├── genome/               # Reference genome files
│   └── metadata/             # Generated metadata files
├── workflow/                  # Workflow components
│   ├── scripts/              # Helper scripts
│   ├── qc_params.snakefile   # QC and parameter generation rules
│   ├── fastp_trimming.snakefile # Trimming rules
│   └── star_alignment.snakefile # Alignment rules
└── Snakefile                  # Main workflow file
```

## Prerequisites

- Snakemake (v6.0+)
- FastQC
- MultiQC
- fastp
- STAR
- samtools

## Setup

1. Clone this repository:
   ```
   git clone https://github.com/yourusername/rna-seq-workflow.git
   cd rna-seq-workflow
   ```

2. Create the necessary directories:
   ```
   mkdir -p resources/config resources/genome resources/metadata
   ```

3. Prepare your sample metadata file (`resources/config/metasheet.csv`):
   ```
   sample,R1,R2
   sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz
   sample2,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz
   ```

4. Configure the workflow parameters in `resources/config/params.yaml`:
   - Adjust memory, threads, and time requirements
   - Set paths to reference genome and annotation files
   - Customize tool-specific parameters

5. Prepare your reference genome:
   - Place your reference genome FASTA in `resources/genome/`
   - Place your GTF annotation file in `resources/genome/`
   - Create a STAR index or use the provided rule:
     ```
     snakemake --cores 16 create_star_index
     ```

## Usage

### Running the Complete Workflow

To run the entire workflow:

```
snakemake --cores <N> --use-conda
```

Replace `<N>` with the number of cores to use.

### Running Specific Stages

The workflow supports running individual stages:

```
# Run just the QC stage
snakemake --cores <N> qc_stage

# Run just the trimming stage
snakemake --cores <N> trimming_stage

# Run just the alignment stage
snakemake --cores <N> alignment_stage
```

### Cluster Execution

For HPC environments, you can use the provided cluster configuration:

```
snakemake --profile slurm
```

Or with explicit cluster parameters:

```
snakemake --cluster "sbatch --mem={resources.mem_mb} --time={params.time} --cpus-per-task={threads}" --jobs 100
```

## Key Features

### Automatic Parameter Generation

The workflow automatically analyzes FastQC results to determine optimal trimming parameters for each sample, stored in `Analysis/QC/Trimming/trimming_params.json`.

### Trimmed Samples Metadata

After trimming, the workflow automatically generates a metadata file (`resources/metadata/trimmed_samples.csv`) with paths to trimmed FASTQ files, which is used for subsequent alignment steps.

### Checkpoint System

The workflow uses checkpoint rules to ensure proper stage completion:
- `qc_complete`: Ensures QC is completed before trimming
- `trimming_complete`: Ensures trimming is completed before alignment
- `alignment_complete`: Marks successful alignment completion

## Output Files

### Quality Control
- `Analysis/QC/FastQC/{sample}/{sample}_R[1,2]_fastqc.html`: FastQC reports
- `Analysis/QC/MultiQC/multiqc_report.html`: MultiQC summary of FastQC results
- `Analysis/QC/Trimming/trimming_params.json`: Sample-specific trimming parameters

### Trimming
- `Analysis/Trimmed/{sample}/{sample}_R[1,2]_trimmed.fastq.gz`: Trimmed FASTQ files
- `Analysis/QC/Trimming/Reports/{sample}_fastp.html`: fastp reports
- `Analysis/QC/Trimming/MultiQC/multiqc_report.html`: MultiQC summary of trimming
- `resources/metadata/trimmed_samples.csv`: Metadata file with trimmed FASTQ paths

### Alignment
- `Analysis/Alignment/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam`: Aligned BAM files
- `Analysis/Alignment/STAR/{sample}/{sample}.Aligned.toTranscriptome.out.bam`: Transcriptome-aligned BAM
- `Analysis/Alignment/STAR/{sample}/{sample}.ReadsPerGene.out.tab`: Gene counts
- `Analysis/Alignment/STAR/{sample}/{sample}.SJ.out.tab`: Splice junctions
- `Analysis/Alignment/STAR/{sample}/{sample}.flagstat.txt`: Alignment statistics
- `Analysis/Alignment/MultiQC/multiqc_report.html`: MultiQC summary of alignment

## Troubleshooting

### Common Issues

1. **Missing input files**: Ensure all input files specified in the metasheet exist
2. **Memory issues**: Adjust memory parameters in `params.yaml` for your system
3. **STAR index errors**: Verify the STAR index was created correctly

### Log Files

All log files are stored in the `logs/` directory, organized by tool:
- `logs/fastqc/`: FastQC logs
- `logs/trimming/`: Trimming logs
- `logs/star/`: Alignment logs
- `logs/workflow/`: Workflow stage logs

## Additional Documentation

- [QC README](docs/README_QC.md): Details on the quality control stage
- [Trimming README](docs/README_TRIMMING.md): Details on the trimming stage
- [Alignment README](docs/README_ALIGNMENT.md): Details on the alignment stage

## License

This workflow is available under the MIT License.

## Citation

If you use this workflow in your research, please cite:

```
Author, A. (Year). RNA-seq Analysis Workflow. GitHub Repository. https://github.com/yourusername/rna-seq-workflow
```

## Contact

For questions or issues, please open an issue on GitHub or contact [your email]. 