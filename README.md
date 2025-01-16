# RNA-seq Data Download and Organization

## Overview
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

# Install required package
conda install -c bioconda sra-tools
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

## Next Steps
After successful download:
1. Proceed to quality control analysis
2. Continue with read trimming
3. Move to alignment and quantification

## Notes
- Keep SRR numbers in a separate file for better project documentation
- Consider backing up the metasheet
- Document any modifications to the standard structure
- Test with a small subset of samples first

## Next Tasks

 - [ ] Add md5 validator for each sample to `download_report.txt`
 - [ ] Add the .yml file for environment
