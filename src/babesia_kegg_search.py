import pandas as pd
from Bio.KEGG import REST
import io
import os
import argparse

def main():
    parser = argparse.ArgumentParser(description='Extract KEGG genes for organisms')
    parser.add_argument('--search_term', '-s', default='babesia', help='Organism search term (default: babesia)')
    parser.add_argument('--output_dir', '-o', default='../data', help='Output directory (default: ../data)')
    args = parser.parse_args()

    # Get all organisms from KEGG
    result = REST.kegg_list("organism").read()
    df = pd.read_table(io.StringIO(result), header=None, sep="\t")

    # Filter by given search term
    filtro = df[df[2].str.contains(args.search_term, case=False, na=False)]
    organism_codes = filtro[1].tolist()

    # Extract genes for each organism (only first 5)
    all_genes_data = []

    for org_code in organism_codes:
        try:
            # Get genes and extract just the IDs
            result = REST.kegg_list(org_code).read()
            genes_df = pd.read_table(io.StringIO(result), header=None, sep="\t")
            
            # Extract first 5 KEGG IDs
            kegg_ids = genes_df[0].head(5).tolist()
            
            # Add each row
            for kegg_id in kegg_ids:
                all_genes_data.append([org_code, kegg_id])
                
        except Exception as e:
            print(f"Skipping {org_code}: {e}")
            continue

    # Create output directory and save single tsv file
    os.makedirs(args.output_dir, exist_ok=True)

    # Save single tsv with 3-letter codes in first column
    if all_genes_data:
        output_file = os.path.join(args.output_dir, f"{args.search_term}_genes.tsv")
        pd.DataFrame(all_genes_data).to_csv(output_file, sep='\t', header=False, index=False)
        print(f"Saved {len(all_genes_data)} genes to {output_file}")

if __name__ == "__main__":
    main()