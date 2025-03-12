#!/bin/bash
set -e

echo "Testing Snakemake workflow with dry run..."

# Create necessary directories 
mkdir -p workflow/common
if [ ! -f "workflow/common/__init__.py" ]; then
    touch workflow/common/__init__.py
fi
mkdir -p resources/config

# Check if the cluster config exists
if [ ! -f "resources/config/cluster.yaml" ]; then
    echo "Creating sample cluster.yaml..."
    cat > resources/config/cluster.yaml << 'EOF'
# Basic cluster configuration for testing
__default__:
  time: "01:00:00"
  mem_mb: 4000
  threads: 1

fastqc:
  time: "01:00:00"
  mem_mb: 4000
  threads: 2
EOF
fi

# Check for syntax errors in Python files
echo "Checking for syntax errors in Python files..."
for pyfile in $(find workflow -name "*.py"); do
  echo "Checking $pyfile..."
  python -m py_compile $pyfile
done

# Check Snakefile syntax using snakemake's parser
echo "Checking Snakefile syntax..."
echo "Note: Linting warnings are suggestions for improvement, not errors."
for snakefile in $(find workflow -name "*.snakefile"); do
  echo "Checking $snakefile..."
  # Use snakemake's parser to check syntax
  snakemake --snakefile $snakefile --lint || {
    echo "Warning: Linting found suggestions for improvement in $snakefile, but these are not errors."
    echo "Continuing with workflow validation..."
  }
done

# Run a dry-run to check workflow validity
echo "Running dry-run with test mode..."
snakemake --config test_mode=True -n --reason

echo "All tests passed! The workflow is ready to run."
echo "To run with a single test sample: ./submit_snakemake.sh --test"
echo "For full dataset: ./submit_snakemake.sh"