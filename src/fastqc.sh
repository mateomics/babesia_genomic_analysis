#!/usr/bin/env bash 
set -e # Only to ensure scrpt executions
set -u # To avoid undefined variables usage
set -o pipefail # To avoid failed runs

#Make the fastqc reports for fastq files 
#Arguments: 
#   $1: directory path with the SRR directories files 
#   $2: Optional argument for a specific output dir   
#The script generate one directory with all the fastqc report on zip and html format

# Check that a path argument was provided
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <path_to_fastq_data>"
    exit 1
fi

data_path="$1"


#if the user specified a output dir 
if [[ $# > 2 ]]; then 
    output_dir=$2 
else 
    output_dir="$data_path"/../../results/fastqc
fi 

# Create a directory for FastQC results"
mkdir -p "$output_dir"

# Iterate over all SRR directories inside the data path
for dir in "$data_path"/*/; do
    # Run FastQC on all FASTQ files in the directory
    fastqc -t 4 -o "$output_dir" "$dir"/*fastq*
done
