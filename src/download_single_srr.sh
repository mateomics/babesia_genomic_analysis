#!/bin/bash

# Ensure robustness
set -e # Only to ensure scrpt executions
set -u # To avoid undefined variables usage
set -o pipefail # To avoid failed runs

output_dir="$1" # Specified on Download.py's 'prefetch_and_conversion()' function
srr_id="$2" # Single SRR ID to download

# Download
prefetch --output-directory "$output_dir" "$srr_id"

srr_dir="$output_dir"/"$srr_id"
# FastQ conversion
fastq-dump "$srr_id" -O "$srr_dir" --split-files --gzip # Separate forward and reverse files, and compress the output

# Iterate over files in the output directory for this SRR
for file in "$srr_dir"/*.fastq.gz; do
    # 1. Verify file exists
    if [ ! -s "$file" ]; then # True if file does not exists or has 0 bytes in it
        echo "$srr_id has been deleted, since it didn't have any content" >> output/summary.txt
        rm "$file"
        continue 1
    fi

    # 2. Verify file has 4 lines -> Identifier, sequence, separator and quality Phred score
    file_lines=$(zcat "$file" | wc -l | cut -d' ' -f1)
    if ((file_lines % 4 != 0)); then # (()) to indicate aritmetic operations
        echo "$srr_id has been deleted, since it was probably truncated" >> output/summary.txt
        rm "$file"
        continue 1
    fi

    # 3. Verify file has as much of '@'s as header lines
    n_lines=$((file_lines / 4))
    line_headers=$(zcat "$file" | grep -c '^@')
    if [ "$n_lines" -ne "$line_headers" ]; then
        echo "$srr_id has been deleted, since it did not have the expected format" >> output/summary.txt
        rm "$file"
        continue 1
    fi
done

