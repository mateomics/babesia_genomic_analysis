#!/usr/bin/env bash 
: '
Make the fastqc reports for fastq files 
Arguments: 
    $1 — directory path with the SRR directories files 

The script generate one directory with all the fastqc report on zip and html format
'
set -e
set -u
set -o pipefail 

# Check that a path argument was provided
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <path_to_fastq_data>"
    exit 1
fi

data_path="$1"

# Create a directory for FastQC results
output_dir="${data_path}/fastqc_results"
mkdir -p "$output_dir"

# Iterate over all SRR directories inside the data path
for dir in "$data_path"/*/; do
    # Run FastQC on all FASTQ files in the directory
    fastqc -t 4 -o "$output_dir" "$dir"/*fastq*
done