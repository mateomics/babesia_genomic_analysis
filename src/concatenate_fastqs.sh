#!/bin/bash
set -e # Only to ensure scrpt executions
set -u # To avoid undefined variables usage
set -o pipefail # To avoid failed runs


output_file="$1"
shift # Once the output file is saved, remove first argument by displacing it
input_files=("$@") # Remaining arguments are now only input files

# Concatenate and re-compress all argument files, one by one
zcat "${input_files[@]}" | gzip > "$output_file" # Instead of using 'subprocess.Popen()' function