"""
Functional enrichment utils module using KEGG
Specialized for B. bovis gene analysis
"""
import io
import os

from typing import Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from Bio.KEGG import REST


def write(text: str, file: str = "enrichment.log", out_dir: str = "../results/EN") -> None:
    """
    Function to write text onto a specified file, log file is set default
    Parameters
    ----------
    -text: str
        The string to be writen 
    -file: str 
        Output file (default "enrichment.log")
    -out_dir: str 
        Output file's directory (default "../results/enrichment")
        
    Returns
    -------
    None
    
    Raises
    ------    
    """ 
    # Make sure path already exist
    dir = os.path.join(out_dir)
    if not os.path.exists(dir):
        os.mkdir(dir) 
    
    # Get file's path 
    summaryfile=os.path.join(dir, file)
    
    # Write text 
    with open (summaryfile, "a") as sum: # Append/add content
        sum.write(f"{text}\n")

def convert_ncbi_to_kegg(ncbi_ids: List[str], organism_code: str = "bbo") -> pd.DataFrame:
    """
    Gets KEGG equivalent ID from an NCBI gene
    
    Parameters
    ----------
    -ncbi_ids: List[str]
        Lista de IDs de NCBI Gene
    -organism_code: str
        Organism's pattern KEGG (default "bbv")
    
    Returns
    -------
    -pd.DataFrame
        Table with NCBI and KEGG's IDs
    """
    try:
        # From NCBI to KEGG
        kegg_string = "+".join([f"{organism_code}:{gene_id}" for gene_id in ncbi_ids]) # "bbv:<n1>+bbv:<n2>+..."
        result = REST.kegg_conv("ncbi-geneid", kegg_string)
        
        content = result.read()
        if not content.strip(): # If did not get anything
            return pd.DataFrame(columns=['ncbi_id', 'kegg_id'])
        
        df = pd.read_table(io.StringIO(content), header=None, sep="\t")
        df.columns = ['kegg_id', 'ncbi_id']  # Columns order inverted to match syntax

        # Pruning IDs, so only basenames are kept
        df['ncbi_id'] = df['ncbi_id'].str.replace('ncbi-geneid:', '')
        df['kegg_id'] = df['kegg_id'].str.replace(f'{organism_code}:', '')
        return df
        
    except Exception as e:
        write(f"Error. Could not get NCBI's IDs from given KEGG's list: {e}")
        return pd.DataFrame(columns=['ncbi_id', 'kegg_id'])

def convert_kegg_ids(id_list: List[str], target_database: str = "uniprot") -> pd.DataFrame:
    """
    Converts KEGG IDs to other database identifiers
    
    Parameters
    ----------
    -id_list: List[str]
        List of KEGG IDs
    -target_database: str
        Target database ("uniprot", "ncbi-geneid", etc.)
    
    Returns
    -------
    -pd.DataFrame
        ID conversion table
        
    Raises
    ------
    Exception
        If ID conversion fails
    """
    try:
        ids_string = "+".join(id_list)
        result = REST.kegg_conv(target_database, ids_string)
        df = pd.read_table(io.StringIO(result.read()), header=None, sep="\t")
        df.columns = ['kegg_id', f'id_{target_database}']
        return df
    except Exception as e:
        write(f"Error in ID conversion: {e}")
        return pd.DataFrame()

def get_kegg_organisms() -> pd.DataFrame:
    """
    Retrieves all available organisms from KEGG database
    
    Returns
    -------
    -pd.DataFrame
        Table with organism codes and descriptions
        
    Raises
    ------
    Exception
        If KEGG API request fails
    """
    try:
        result = REST.kegg_list("organism")
        df = pd.read_table(io.StringIO(result.read()), header=None, sep="\t")
        df.columns = ['code', 'description']
        return df
    except Exception as e:
        write(f"Error retrieving organisms: {e}")
        return pd.DataFrame()

def search_organism(search_term: str = "bovis") -> pd.DataFrame:
    """
    Finds organisms matching the search term
    
    Parameters
    ----------
    -search_term : str
        Search term (default "bovis")
    
    Returns
    -------
    -pd.DataFrame
        Organisms matching the search criteria
        
    Raises
    ------
    Exception
        If search operation fails
    """
    try:
        organisms = get_kegg_organisms()
        if organisms.empty:
            write("Organisms absent on KEGG")
            return pd.DataFrame()
        
        # Only organisms that matches the given pattern
        mask = organisms['description'].str.contains(search_term, case=False, na=False) # Case insensitive, Na's as false
        return organisms[mask]
    except Exception as e:
        write(f"Error in organism search: {e}")
        return pd.DataFrame()

def get_organism_genes(organism_code: str = "bbo") -> pd.DataFrame:
    """
    Retrieves all genes for a specific organism
    
    Parameters
    ----------
    -organism_code: str
        Organism code (default "bbv" for B. bovis)
    
    Returns
    -------
    -pd.DataFrame
        Genes from the specified organism
        
    Raises
    ------
    Exception
        If gene retrieval fails
    """
    try:
        result = REST.kegg_list(organism_code)
        df = pd.read_table(io.StringIO(result.read()), header=None, sep="\t")
        df.columns = ['gene_id', 'description']
        return df
    except Exception as e:
        write(f"Error retrieving genes for {organism_code}: {e}")
        return pd.DataFrame()

def get_organism_pathways(organism_code: str = "bbo") -> pd.DataFrame:
    """
    Retrieves metabolic pathways for a specified organism
    
    Parameters
    ----------
    -organism_code: str
        Organism cod (default "bbv")
    
    Returns
    -------
    -pd.DataFrame
        Metabolic pathways for the organism
        
    Raises
    ------
    Exception
        If pathway retrieval fails
    """
    try:
        result = REST.kegg_list("pathway", organism_code)
        df = pd.read_table(io.StringIO(result.read()), header=None, sep="\t")
        df.columns = ['pathway_id', 'description']
        return df
    
    except Exception as e:
        write(f"Error retrieving pathways for {organism_code}: {e}")
        return pd.DataFrame()

def get_gene_info(gene_id: str) -> str:
    """
    Retrieves detailed information for a specified gene, using KEGG_get()
    
    Parameters
    ----------
    -gene_id: str
        Gene ID (for example, "bbv:BBOV_III000110")
    
    Returns
    -------
    -str
        Complete gene information
        
    Raises
    ------
    Exception
        If gene information retrieval fails
    """
    try:
        result = REST.kegg_get(gene_id)
        return result.read()
    except Exception as e:
        write(f"Error retrieving information for {gene_id}: {e}")
        return ""

def search_genes_by_term(search_term: str, organism_code: str = "genes") -> pd.DataFrame:
    """
    Finds genes matching the search term
    
    Parameters
    ----------
    -search_term : str
        Search term
    -organism_code : str, optional
        Organism code for filtering
    
    Returns
    -------
    -pd.DataFrame
        Genes matching the search criteria
        
    Raises
    ------
    Exception
        If gene search fails
    """
    try:
        result = REST.kegg_find(organism_code, search_term)
        df = pd.read_table(io.StringIO(result.read()), header=None, sep="\t")
        df.columns = ['gene_id', 'description']
        return df
    except Exception as e:
        write(f"Error in gene search: {e}")
        return pd.DataFrame()

def get_pathways_by_gene(gene_id: str) -> pd.DataFrame:
    """
    Retrieves metabolic pathways associated with a gene
    
    Parameters
    ----------
    -gene_id: str
        Gene ID
    
    Returns
    -------
    -pd.DataFrame
        Pathways associated with the gene
        
    Raises
    ------
    Exception
        If pathway association retrieval fails
    """
    try:
        result = REST.kegg_link("pathway", gene_id)
        df = pd.read_table(io.StringIO(result.read()), header=None, sep="\t")
        df.columns = ['gene_id', 'pathway_id']
        return df
    except Exception as e:
        write(f"Error retrieving pathways for {gene_id}: {e}")
        return pd.DataFrame()

def get_genes_by_pathway(pathway_id: str) -> pd.DataFrame:
    """
    Retrieves genes associated with a metabolic pathway
    
    Parameters
    ----------
    -pathway_id: str
        Pathway ID (for example, "path:bbv00010")
    
    Returns
    -------
    -pd.DataFrame
        Genes associated with the pathway
        
    Raises
    ------
    Exception
        If gene association retrieval fails
    """
    try:
        organism_code = pathway_id.replace("path:", "").split("0")[0]
        result = REST.kegg_link(organism_code, f"path:{pathway_id}")
        df = pd.read_table(io.StringIO(result.read()), header=None, sep="\t")
        df.columns = ['pathway_id', 'gene_id']
        return df
    except Exception as e:
        write(f"Error retrieving genes for {pathway_id}: {e}")
        return pd.DataFrame()

def get_pathway_image(pathway_id: str, organism_code: str = "map") -> bytes:
    """
    Retrieves pathway image for visualization
    
    Parameters
    ----------
    -pathway_id : str
        Pathway ID (e.g., "00061" for fatty acid biosynthesis)
    -organism_code : str
        Organism code for organism-specific pathways (default "map")
    
    Returns
    -------
    -bytes
        Image data in bytes
        
    Raises
    ------
    Exception
        If image retrieval fails
    """
    try:            
        result = REST.kegg_get(f"{organism_code}{pathway_id}", "image")
        return result.read()
    except Exception as e:
        write(f"Error retrieving pathway image for {pathway_id}: {e}")
        return b""

def get_kegg_databases() -> pd.DataFrame:
    """
    Retrieves list of all available KEGG databases
    
    Returns
    -------
    -pd.DataFrame
        Table with database codes and descriptions
        
    Raises
    ------
    Exception
        If database retrieval fails
    """
    try:
        result = REST.kegg_list("database")
        df = pd.read_table(io.StringIO(result.read()), header=None, sep="\t")
        df.columns = ['database_code', 'description']
        return df
    except Exception as e:
        write(f"Error retrieving KEGG databases: {e}")
        return pd.DataFrame()

def get_database_info(database_code: str) -> str:
    """
    Retrieves information about a specific KEGG database
    
    Parameters
    ----------
    -database_code: str
        Database code (for example: "pathway", "compound")
    
    Returns
    -------
    -str
        Database information and statistics
        
    Raises
    ------
    Exception
        If database info retrieval fails
    """
    try:
        result = REST.kegg_info(database_code)
        return result.read()
    except Exception as e:
        write(f"Error retrieving info for database {database_code}: {e}")
        return ""

def perform_enrichment_analysis(gene_list: List[str], organism_code: str = "bbo") -> pd.DataFrame:
    """
    Performs functional enrichment analysis for a gene list
    
    Parameters
    ----------
    -gene_list: List[str]
        List of gene IDs
    -organism_code: str
        Organism code
    
    Returns
    -------
    -pd.DataFrame
        Enrichment analysis results by pathway
        
    Raises
    ------
    Exception
        If enrichment analysis fails
    """
    try:

        gene_list_with_prefix = [f"{organism_code}:{gene}" for gene in gene_list]

        organism_pathways = get_organism_pathways(organism_code)
        
        enrichment_results = [] # Empty list

        for _, pathway in organism_pathways.iterrows():
            pathway_id = pathway['pathway_id']
            pathway_genes = get_genes_by_pathway(pathway_id)
            
            if pathway_genes.empty:
                continue
            genes_in_pathway = set(pathway_genes['gene_id']) # In both cases to avoid duplicates
            genes_of_interest = set(gene_list_with_prefix) # Use the prefixed list
            
            # Only desired genes also present at pathway
            common_genes = [present for present in genes_in_pathway if present in genes_of_interest]
                
            if len(common_genes) == 0:
                continue

            # Creating results list
            enrichment_results.append({
                'pathway_id': pathway_id,
                'pathway_description': pathway['description'],
                'total_genes_in_pathway': len(genes_in_pathway),
                'genes_of_interest_in_pathway': len(common_genes),
                'percentage_in_pathway': (len(common_genes) / len(genes_of_interest)) * 100,
                'common_genes': ', '.join(common_genes)
            })
        
        (df_enrichment := pd.DataFrame(enrichment_results)).to_csv("../results/EN/enrichment_results.csv", sep="\t", index=False)
        return df_enrichment
    
    except Exception as e:
        write(f"Error in enrichment analysis: {e}")
        return pd.DataFrame()

def process_gene_table(gene_table: pd.DataFrame) -> List[str]:
    """
    Processes a gene table for enrichment analysis
    
    Parameters
    ----------
    -gene_table: pd.DataFrame
        Table with gene information
    
    Returns
    -------
    -List[str]
        List of unique gene IDs
        
    Raises
    ------
    Exception
        If table processing fails
    """
    try:   
        unique_genes = gene_table.iloc[:, -1].dropna().unique().tolist()
        return [gene for gene in unique_genes if gene and str(gene) != 'nan']
    
    except Exception as e:
        print(f"Error processing table: {e}")
        return []

def plot_enrichment_bars(enrichment_df: pd.DataFrame,
                        title: str = "Functional Enrichment Analysis",
                        output_file: str = "enrichment_barplot.png",
                        out_dir: str = "../results/EN",
                        max_pathways: int = 15) -> str:
    """
    Generates a clean barplot for enrichment analysis results
    
    Parameters
    ----------
    -enrichment_df: pd.DataFrame
        DataFrame with enrichment results
    -title: str
        Plot title
    -output_file: str
        Output filename
    -out_dir: str
        Output directory
    -max_pathways: int
        Maximum number of pathways to display
    
    Returns
    -------
    -str
        Path to saved file
    """
    try:
        if enrichment_df.empty:
            return ""
        
        os.makedirs(out_dir, exist_ok=True)
        
        # Sort and select top pathways
        plot_data = enrichment_df.sort_values('percentage_in_pathway', ascending=True)
        if len(plot_data) > max_pathways:
            plot_data = plot_data.tail(max_pathways)
        
        # Create shorter pathway descriptions with 3-letter organism codes
        short_descriptions = []
        for desc in plot_data['pathway_description']:
            # Extract organism and convert to 3-letter code
            parts = desc.split(' - ')
            pathway_name = parts[0]
            
            # Get organism code (last part) and convert to 3 letters
            if len(parts) > 1:
                organism = parts[-1]
                # Take first 3 letters of organism and capitalize
                org_code = organism[:3].upper()
                short_desc = f"{pathway_name} ({org_code})"
            else:
                short_desc = pathway_name
            
            # Truncate if still too long
            if len(short_desc) > 60:
                short_desc = short_desc[:57] + '...'
            short_descriptions.append(short_desc)
        
        # Create plot
        plt.figure(figsize=(12, 8)) # Increased figure size for better appreciation
        
        y_pos = range(len(plot_data))
        
        # Increase spacing between bars
        bars = plt.barh(y_pos, plot_data['percentage_in_pathway'], 
                       color='darkblue', alpha=0.8, height=0.6) # Reduced height for more spacing
        
        # Format labels with increased spacing
        plt.yticks(y_pos, short_descriptions, fontsize=10)
        plt.xlabel('Gene Percentage (%)', fontsize=12, fontweight='bold')
        plt.ylabel('Pathways', fontsize=12, fontweight='bold')
        plt.title(title, fontsize=14, pad=20, fontweight='bold')
        
        # Add value labels on bars
        for i, (bar, value) in enumerate(zip(bars, plot_data['percentage_in_pathway'])):
            plt.text(bar.get_width() + 0.3, bar.get_y() + bar.get_height()/2, 
                    f'{value:.1f}%', ha='left', va='center', fontsize=9, fontweight='bold')
        
        # Improve layout with better spacing
        plt.tight_layout()
        plt.grid(axis='x', alpha=0.3, linestyle='--')
        
        # Adjust subplot parameters to make plot more prominent
        plt.subplots_adjust(left=0.3, right=0.9, top=0.9, bottom=0.1)
        
        # Save with high quality
        output_path = os.path.join(out_dir, output_file)
        plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='lightgray')
        plt.close()
        
        write(f"Plot saved: {output_path}")
        return output_path
        
    except Exception as e:
        write(f"Error generating barplot: {e}")
        return ""