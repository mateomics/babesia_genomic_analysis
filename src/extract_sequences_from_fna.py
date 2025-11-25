#!/usr/bin/env python3
"""
Extracts the secuence of the genome given coordinates of a TSV file

Usage:
    python extract_sequences_simple.py -g <genome_file> -c <coords_file> -o <output>

"""

from Bio import SeqIO
import pandas as pd
import argparse
import os


def validate_files(genome_file: str, coords_file: str) -> None:
    """
    Validate that input files exist and are not empty.
    
    Args:
        genome_file (str): Path to the genome FASTA file
        coords_file (str): Path to the coordinates TSV file
    
    Raises:
        FileNotFoundError: If files don't exist
        ValueError: If files are empty
    """
    # Check genome file
    if not os.path.exists(genome_file):
        raise FileNotFoundError(f"Genome file not found: {genome_file}")
    if os.path.getsize(genome_file) == 0:
        raise ValueError(f"Genome file is empty: {genome_file}")
    
    # Check coords file
    if not os.path.exists(coords_file):
        raise FileNotFoundError(f"Coordinates file not found: {coords_file}")
    if os.path.getsize(coords_file) == 0:
        raise ValueError(f"Coordinates file is empty: {coords_file}")
    
    print("[OK] Input files validated successfully")


def validate_coords_dataframe(coords_df: pd.DataFrame) -> None:
    """
    Validate that the coordinates dataframe has required columns and valid data.
    
    Args:
        coords_df (pd.DataFrame): DataFrame containing genomic coordinates
    
    Raises:
        ValueError: If required columns are missing or data is invalid
    """
    required_columns = ['geneid', 'chr', 'start', 'end', 'strand']
    
    # Check for required columns
    missing_columns = [col for col in required_columns if col not in coords_df.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns: {', '.join(missing_columns)}")
    
    # Check for empty dataframe
    if len(coords_df) == 0:
        raise ValueError("Coordinates file contains no data rows")
    
    # Validate start < end for all rows
    invalid_coords = coords_df[coords_df['start'] >= coords_df['end']]
    if len(invalid_coords) > 0:
        print(f"[WARNING] Found {len(invalid_coords)} genes with start >= end:")
        for idx, row in invalid_coords.iterrows():
            print(f"  - {row['geneid']}: start={row['start']}, end={row['end']}")
        raise ValueError("Invalid coordinates detected: start must be less than end")
    
    print("[OK] Coordinates dataframe validated successfully")


def coords_to_fasta(genome_file:str, coords_file:str, output_file:str) -> None:
    """
    Extract sequences from a genome file based on coordinates from a TSV file.
    
    Reads a genome FASTA file and a TSV file containing genomic coordinates,
    then extracts the sequences at those coordinates and writes them to an output
    FASTA file. Handles reverse complement for genes on the negative strand.
    
    Args:
        genome_file (str): Path to the genome FASTA file
        coords_file (str): Path to the TSV file with columns: geneid, chr, start, end, strand, DE_status
        output_file (str): Path where the output FASTA file will be written
    
    Returns:
        None: Writes results to the output file and prints summary statistics
    """

    # Validate input files
    validate_files(genome_file, coords_file)
    
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir)
            print(f"[INFO] Created output directory: {output_dir}")
        except OSError as error:
            raise IOError(f"Error creating output directory: {error}")
 
    print("[INFO] Loading genome...")
    genome = {}
    try:
        for record in SeqIO.parse(genome_file, "fasta"):
            genome[record.id] = record.seq
            print(f"  - {record.id}: {len(record.seq):,} bp")
    except Exception as error:
        raise IOError(f"Error reading genome file: {error}")
    
    # Check if genome was loaded successfully
    if not genome:
        raise ValueError("No sequences found in genome file")
        
    print("\n[INFO] Loading coordinates...")
    try:
        coords_df = pd.read_csv(coords_file, sep='\t')
    except Exception as error:
        raise IOError(f"Error reading coordinates file: {error}")
    
    # Validate coordinates dataframe
    validate_coords_dataframe(coords_df)
    
    print(f"Total genes to extract: {len(coords_df)}")
    
    print("\nGathering sequences...")
    extracted = 0
    errors = 0
    
    try:
        with open(output_file, 'w') as output_fasta:
            for idx, row in coords_df.iterrows():
                geneid = row['geneid']
                chromosome = row['chr']
                start = int(row['start'])
                end = int(row['end'])
                strand = row['strand']
                de_status = row['DE_status'] if 'DE_status' in row else 'N/A'
                
                # Verify that the chromosome exists
                if chromosome not in genome:
                    print(f"[ERROR] Chromosome not found: {chromosome}")
                    errors += 1
                    continue
                
                # Extract the sequence
                chrom_seq = genome[chromosome]
                start_idx = start - 1  # Convert to 0-based indexing
                end_idx = end
                
                sequence = chrom_seq[start_idx:end_idx]
                
                # Gets the reverse complement if the strand is negative
                if strand == '-':
                    sequence = sequence.reverse_complement()
                
                # Creates the header of the file
                header = f">{geneid}|{chromosome}|{start}-{end}|strand_{strand}|{de_status}"
                
                # Writes the header to the output file
                output_fasta.write(f"{header}\n")
                
                # Writes the sequence with 80 chars per line following fasta standards
                seq_str = str(sequence)
                for i in range(0, len(seq_str), 80):
                    output_fasta.write(seq_str[i:i+80] + '\n')
                
                extracted += 1
                
                if extracted % 10 == 0:
                    print(f"[INFO] Extracted: {extracted}")
    except IOError as error:
        raise IOError(f"Error writing to output file: {error}")
    
    print(f"\n{'='*50}")
    print(f"[INFO] SUMMARY:")
    print(f"[INFO] Extracted sequences: {extracted}")
    print(f"[INFO] Errors: {errors}")
    print(f"[INFO] Output file path: {output_file}")
    print(f"[OK]")
    print(f"{'='*50}")

def parser() -> argparse.ArgumentParser: 
    """
    Function to parser the arguments provided by the user 
    Args: 
    Returns: 
        -args(): object that store all the arguments parsed 
    """  

    # Initialize the object parser 
    parser= argparse.ArgumentParser(description="Parser to store the user information") 

    #arguments to be specified
    parser.add_argument("--genome-file", "-g",
                        default=None,
                        type=str,
                        help="Path of the genome fasta file")
    parser.add_argument("--coords-file", "-c",
                        default=None,
                        type=str,
                        help="Path of the coordinates file")
    parser.add_argument("--output-fasta", "-o",
                        default="./results/output_fasta.fasta",
                        type=str,
                        help="Path of the output fasta file")
    
    # Get the arguments parsed 
    args= parser.parse_args() 
    
    return args

def main():

    arguments = parser()

    genome_file = arguments.genome_file
    coords_file = arguments.coords_file
    output_file = arguments.output_fasta

    try:
        coords_to_fasta(genome_file, coords_file, output_file)
    except (FileNotFoundError, ValueError, IOError) as error:
        print(f"\n[ERROR] {error}")
        exit(1)

if __name__ == "__main__":
    main()
