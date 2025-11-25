#!/usr/bin/env python3
"""
Script to add a GeneID column to the DE_genes_in_bovis.tab file
based on matches with the Babesia_bovis.old file
"""

import re
import subprocess
from pathlib import Path
from typing import Dict
import argparse

def parse_fasta_headers(fasta_file: Path) -> Dict[str, str]:
    """
    Reads the FASTA file and creates a dictionary mapping ID -> GeneID
    """
    id_to_geneid = {}
    
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                # Extract the ID (after >lcl| and before the space)
                match_id = re.search(r'>lcl\|([^\s]+)', line)
                # Extract the GeneID
                match_geneid = re.search(r'GeneID:(\d+)', line)
                
                if match_id and match_geneid:
                    sequence_id = match_id.group(1)
                    gene_id = match_geneid.group(1)
                    id_to_geneid[sequence_id] = gene_id
    
    return id_to_geneid

def add_geneid_column(tab_file: Path, output_file: Path, id_to_geneid: Dict[str, str]) -> None:
    """
    Reads the .tab file, adds the GeneID column, and saves the result
    """
    with open(tab_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            line = line.rstrip('\n')
            fields = line.split('\t')
            
            if len(fields) >= 2:
                # The second column (index 1) contains the ID
                sequence_id = fields[1]
                
                # Search for the corresponding GeneID
                gene_id = id_to_geneid.get(sequence_id, 'NA')
                
                # Add the GeneID column to the end
                new_line = line + '\t' + gene_id
                outfile.write(new_line + '\n')
            else:
                # If the line doesn't have enough columns, write it as is
                outfile.write(line + '\n')

def run_blastn(query_file: Path, db_path: str, output_file: Path) -> None:
    """
    Executes BLASTN using subprocess and waits for completion
    
    Args:
        query_file: Path to the query FASTA file
        db_path: Path to the BLAST database
        output_file: Path where BLAST results will be saved
    """
    # Create output directory if it doesn't exist
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    # BLASTN command
    blastn_cmd = [
        "blastn",
        "-query", str(query_file),
        "-db", db_path,
        "-out", str(output_file),
        "-outfmt", "6 qseqid sseqid pident gaps length qlen slen qcovs evalue score"  # Tabular format
    ]
    
    print(f"[INFO] Running BLASTN...")
    print(f"[INFO] Command: {' '.join(blastn_cmd)}")
    
    try:
        # Execute BLASTN and wait for completion
        result = subprocess.run(blastn_cmd, check=True, capture_output=True, text=True)
        print(f"[INFO] BLASTN completed successfully")
        if result.stdout:
            print(f"[INFO] STDOUT: {result.stdout}")
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] BLASTN failed with error code {e.returncode}")
        print(f"[ERROR] STDERR: {e.stderr}")
        raise
    except FileNotFoundError:
        print("[ERROR] blastn command not found. Make sure BLAST+ is installed and in PATH")
        raise

def parser() -> argparse.ArgumentParser:
    """
    Function to parse the arguments provided by the user
    
    Returns:
        argparse.Namespace: Object that stores all the parsed arguments
    """
    # Initialize the object parser
    parser = argparse.ArgumentParser(
        description="Add GeneID column to tab file based on FASTA file"
    )
    
    # Arguments to be specified
    parser.add_argument("--fasta-file", "-f",
                        default=None,
                        type=str,
                        help="Path to the FASTA file with GeneID annotations")
    parser.add_argument("--query-file", "-q",
                        default=None,
                        type=str,
                        help="Path to the query FASTA file for BLAST (e.g., DE_genes_sequences.fasta)")
    parser.add_argument("--db-path", "-d",
                        default=None,
                        type=str,
                        help="Path to the BLAST database")
    parser.add_argument("--output-file", "-o",
                        default="./results/output_with_geneid.tab",
                        type=str,
                        help="Path to the output tab file")
    
    # Get the arguments parsed
    args = parser.parse_args()
    
    return args

def main() -> None:
    arguments = parser()
    
    # Define file paths from arguments
    fasta_file = Path(arguments.fasta_file)
    query_file = Path(arguments.query_file)
    db_path = arguments.db_path
    output_file = Path(arguments.output_file)
    
    # Create output directory if it doesn't exist
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    # Define BLAST output path (temporary file in the same directory as output_file)
    blast_output = output_file.parent/"blast_results.tab"
    
    # Step 1: Run BLASTN
    print(f"[INFO] Running BLASTN")
    print(f"[INFO] Query file: {query_file}")
    print(f"[INFO] Database: {db_path}")
    run_blastn(query_file, db_path, blast_output)
    print(f"[OK] BLASTN executed")
    
    # Step 2: Parse FASTA file for GeneID mappings
    print(f"\n[INFO] Reading FASTA file: {fasta_file}")
    id_to_geneid = parse_fasta_headers(fasta_file)
    print(f"[INFO] Found {len(id_to_geneid)} ID -> GeneID mappings")
    print(f"[OK] FASTA file processed")
    
    # Step 3: Add GeneID column to BLAST results
    print(f"\n[INFO] Processing BLAST results: {blast_output}")
    add_geneid_column(blast_output, output_file, id_to_geneid)
    print(f"[OK] Output file created: {output_file}")

    
    # Show some statistics
    matched = 0
    total = 0
    with open(output_file, 'r') as f:
        for line in f:
            total += 1
            if not line.strip().endswith('\tNA'):
                matched += 1
    
    print(f"\n[INFO] Statistics:")
    print(f"[INFO] Total lines: {total}")
    print(f"[INFO] Matches found: {matched}")
    print(f"[INFO] No match: {total - matched}")

if __name__ == "__main__":
    main()
