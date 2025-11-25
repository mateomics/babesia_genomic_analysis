"""
Enrichment analysis pipeline for B. bovis
"""

import argparse
import pandas as pd
import os
from utils_en import (
    process_gene_table,
    convert_ncbi_to_kegg, 
    perform_enrichment_analysis,
    plot_enrichment_bars,
    write)

def main():
    # Define arguments
    parser = argparse.ArgumentParser(description='Functional enrichment analysis for B. bovis genes')
    
    parser.add_argument(
        '--input', 
        '-i', 
        required=True,
        help='Input gene table file (tab as field separator)')
    parser.add_argument(
        '--output_dir', 
        '-o', 
        default='../results/EN',
        help='Output directory for results (default: ../results/EN)')
    parser.add_argument(
        '--organism',
        '-org',
        default='bbo',
        help='KEGG organism code (default "bbo" for B. bovis)')
    parser.add_argument(
        '--kegg_ids',
        '-k',
        action='store_true',
        help='Flag to indicate the input is already a KEGG IDs list')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    write("Starting enrichment analysis pipeline", out_dir=args.output_dir)
    write(f"Input file: {args.input}", out_dir=args.output_dir)
    write(f"Organism: {args.organism}", out_dir=args.output_dir)
    write(f"KEGG IDs mode: {args.kegg_ids}", out_dir=args.output_dir)
    
    try:
        # Load input table
        write("Loading gene table...", out_dir=args.output_dir)
        gene_table = pd.read_csv(args.input, sep='\t', header=None)
        
        write(f"Loaded {len(gene_table)} rows", out_dir=args.output_dir)
        
        # Extract gene IDs
        write("Extracting gene IDs...", out_dir=args.output_dir)
        gene_ids = process_gene_table(gene_table)
        write(f"Found {len(gene_ids)} unique gene IDs", out_dir=args.output_dir)
        
        if not gene_ids:
            write("Error: No gene IDs found", out_dir=args.output_dir)
            return
        
        # Skip conversion if using KEGG IDs directly
        if args.kegg_ids:
            write("Skipping NCBI to KEGG conversion (using direct KEGG IDs)", out_dir=args.output_dir)
            kegg_ids = gene_ids # Use the extracted IDs directly as KEGG IDs
            
            # Just a dummy conversion table for subsequent consistency
            kegg_df = pd.DataFrame({
                'kegg_id': kegg_ids,
                'ncbi_id': ['direct_kegg'] * len(kegg_ids) # Placeholder
            })
            
        else:
            # Normal conversion flow
            write("Converting NCBI to KEGG IDs...", out_dir=args.output_dir)
            kegg_df = convert_ncbi_to_kegg(gene_ids, args.organism)
            write(f"Converted {len(kegg_df)} IDs", out_dir=args.output_dir)
            
            if kegg_df.empty:
                write("No matches: No KEGG ID found", out_dir=args.output_dir)
                return
            
            kegg_ids = kegg_df['kegg_id'].tolist()
        
        # Save conversion table (even if using direct KEGG IDs)
        conversion_file = os.path.join(args.output_dir, "gene_ids_table.csv")
        kegg_df.to_csv(conversion_file, index=False)
        write(f"Gene IDs table saved: {conversion_file}", out_dir=args.output_dir)
        
        # Enrichment analysis
        write("Running enrichment analysis...", out_dir=args.output_dir)
        write(f"Using {len(kegg_ids)} KEGG IDs: {kegg_ids[:5]}...", out_dir=args.output_dir)
        
        enrichment_df = perform_enrichment_analysis(kegg_ids, args.organism)
        write(f"Found {len(enrichment_df)} enriched pathways", out_dir=args.output_dir)
        
        # Save results and generate plots
        if not enrichment_df.empty:
            # Save enrichment results
            results_file = os.path.join(args.output_dir, "enrichment_results.csv")
            enrichment_df.to_csv(results_file, index=False)
            write(f"Results saved: {results_file}", out_dir=args.output_dir)
            
            # Generate plots
            write("Generating plot...", out_dir=args.output_dir)
            plot_enrichment_bars(
                enrichment_df, 
                title=f"Enriched Pathways {args.organism}",
                output_file="enrichment_bars.png",
                out_dir=args.output_dir)
            
        else:
            write("No enriched pathways found", out_dir=args.output_dir)
        
        write("Pipeline completed successfully", out_dir=args.output_dir)
        
    except Exception as e:
        write(f"Pipeline failed: {e}", out_dir=args.output_dir)
        raise

if __name__ == "__main__":
    main()