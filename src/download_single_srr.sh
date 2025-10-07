#!/bin/bash

output_dir="$1" # Specified on Download.py's 'prefetch_and_conversion()' function
srr_id="$2" # Single SRR ID to download

# Download
prefetch --output-directory "$output_dir" "$srr_id"
# FastQ conversion
fastq-dump "$srr_id" -O "$output_dir/$srr_id" --split-files --gzip # Separate forward and reverse files, and compress the output