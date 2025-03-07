# RNA-seq pipeline analysis optimized for Ohio Supercomputer center

## RNA-seq Data Download and Organization

### Overview
This section of the pipeline handles downloading and organizing RNA-seq data from the Sequence Read Archive (SRA). The script is designed to:
- Download samples from SRA using prefetch and fasterq-dump
- Organize samples into lab/project-specific directories
- Generate a metasheet for downstream analysis
- Provide detailed download reports

## Directory Structure
The script will create the following structure:
```
base_directory/
├── Lab1_Sample/
│   ├── SRRxxx_R1.fastq.gz
│   └── SRRxxx_R2.fastq.gz
├── Lab2_Sample/
│   └── [fastq files]
├── resources/
│   └── config/
│       └── metasheet.csv
└── download_report.txt
```

## Setup Instructions

### 1. Environment Setup
```bash
# Load conda module
module load miniconda3

# Create and activate environment
conda create -n your_env_name
conda activate your_env_name
```

### 2. Modify the Download Script

Edit `download_sra.sh` to include your samples:

```bash
# Define your base directory
BASE_DIR="/path/to/your/project"

# Define your samples by group/lab
declare -A lab_samples=(
    ["LabName1"]="SRR1 SRR2"
    ["LabName2"]="SRR1 SRR2"
)
```

Key elements to modify:
- `BASE_DIR`: Set your project directory path
- `lab_samples`: Add your SRR numbers grouped by lab/project
- SLURM parameters: Adjust based on your needs and system

### 3. Run the Download

```bash
sbatch download_sra.sh
```

## Output Files

### 1. Metasheet
The script generates `metasheet.csv` with the structure:
```csv
sample,R1,R2
LabName1_SRR1234567,/path/to/R1.fastq.gz,/path/to/R2.fastq.gz
```
This metasheet is used by downstream analysis steps.

### 2. Download Report
The `download_report.txt` includes:
- Download status for each sample
- File sizes and locations
- Summary statistics
- Error reports (if any)

## Resource Requirements

### Default SLURM Configuration
```bash
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=32G
```
Adjust these based on:
- Number of samples
- Size of fastq files
- System capabilities

## Troubleshooting Guide

### Common Issues and Solutions

1. Download Failures
   - Verify SRR numbers are correct
   - Check internet connectivity
   - Ensure conda environment is activated
   - Verify disk space availability

2. Organization Issues
   - Check write permissions in base directory
   - Verify directory structure exists
   - Check for conflicting files

### Validation Steps
After download completion:
1. Check `download_report.txt`
2. Verify file sizes are reasonable (typically several GB per fastq file)
3. Confirm all expected files are present
4. Validate metasheet entries

## Customization

### Adding New Sample Groups
```bash
# Example addition to lab_samples array
["NewLab"]="SRR9876543 SRR9876544"
```

### Modifying File Organization
The script can be adapted for different organizational needs:
- Change directory naming scheme
- Modify file naming patterns
- Adjust metasheet format

## Notes
- Keep SRR numbers in a separate file for better project documentation
- Consider backing up the metasheet
- Document any modifications to the standard structure
- Test with a small subset of samples first

# Quality Control Analysis
## Overview
This section of the pipeline performs quality control analysis of RNA-seq data using FastQC and MultiQC. The workflow is split into two parts:
- FastQC analysis of individual samples
- MultiQC aggregation of all FastQC reports

## Directory Structure
The quality control workflow creates the following structure:

```
base_directory/
├── Analysis/
│   └── QC/
│       ├── FastQC/
│       │   ├── Lab1_SRRxxx/
│       │   │   ├── SRRxxx_R1_fastqc.html
│       │   │   ├── SRRxxx_R1_fastqc.zip
│       │   │   ├── SRRxxx_R2_fastqc.html
│       │   │   └── SRRxxx_R2_fastqc.zip
│       │   └── Lab2_SRRxxx/
│       │       └── [FastQC files]
│       └── MultiQC/
│           ├── multiqc_report.html
│           └── multiqc_report.pdf
└── logs/
    ├── fastqc/
    │   └── [FastQC log files]
    ├── multiqc/
    │   └── [MultiQC log files]
    └── slurm/
        └── [SLURM output files]
```

## Setup and Execution
### 1. FastQC Analysis
Run FastQC on all samples using Snakemake:
```bash
# Submit FastQC workflow
sbatch submit_fastqc.sh
```

### 2. MultiQC Report Generation
After FastQC completion, generate the MultiQC report:
```bash
# Submit MultiQC workflow
sbatch submit_multiqc.sh
```

## Configuration
### FastQC Parameters
Edit `resources/config/params.yaml` to modify:
```yaml
# FastQC settings
fastqc_threads: 2
fastqc_memory: 4000
fastqc_time: "01:00:00"
```

### SLURM Configuration
Default settings in submit scripts:
```bash
# FastQC job settings
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --mem=4G

# MultiQC job settings
#SBATCH --time=0:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
```

## Output Files
### 1. FastQC Reports
For each sample (R1 and R2):
- HTML report with quality metrics
- ZIP archive containing raw data

### 2. MultiQC Summary
- Interactive HTML report combining all FastQC results
- PDF version for sharing/printing

## Quality Metrics Analyzed
FastQC examines:
- Base sequence quality
- Sequence length distribution
- GC content
- Sequence duplication levels
- Overrepresented sequences
- Adapter content

## Troubleshooting
### Common Issues
1. FastQC Job Failures
   - Check input file paths in metasheet
   - Verify FastQC module is loaded
   - Check resource requirements

2. MultiQC Issues
   - Ensure all FastQC jobs completed
   - Verify FastQC output directory structure
   - Check MultiQC dependencies

### Validation Steps
1. Confirm FastQC completion:
   - Check all expected output files exist
   - Review individual FastQC logs
   
2. Verify MultiQC report:
   - Ensure all samples are included
   - Check for any failed analyses

## Notes
- FastQC and MultiQC workflows are separated for better management
- Increase `--latency-wait` if experiencing file system delays
- Consider adjusting threads and memory based on sample size
```

## Next Tasks

 - [ ] Add md5 validator for each sample to `download_report.txt`
 - [ ] Add the .yml file for environment
 - [ ] 