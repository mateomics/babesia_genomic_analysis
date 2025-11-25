"""
Program to make the differential expression analysis with a count matrix with the gene data in ever experimental condition 


""" 
import pandas as pd 
import argparse
import utils_DE as ut 
import subprocess



def parser()->argparse.ArgumentParser: 
    """
    Function to parse the arguments provided by the user 
    Args: 
    Returns: 
        -args(): object that store all the arguments parsed 
    """  
    #We initialize de object parser 
    parser= argparse.ArgumentParser(description="Parser to store the user information") 

    #arguments to be specified
    parser.add_argument("--matrix",
                        default=None,
                        help="Specifie the path of the count matrix")
    parser.add_argument("--dess",
                        default=None,
                        help="Specifie the path with a design matrix created previously")
    parser.add_argument("--plots",
                        default=3,
                        type=int,
                        help="Specifie the plots do you want separated by a coma plots supported=vulcano,heat map ; 0 for any plot 1 volcano plot, 2 heat map, 3 both plots")   
    parser.add_argument("--pref",
                        default="",
                        help="Specifie a prefix to save the plots generated") 
    parser.add_argument("--logF",
                        default=1,
                        type=float,
                        help="Specifie the tresshold for Diferential expressed genes log2 fold change")
    parser.add_argument("--padj",
                        default=0.05,
                        type=float,
                        help="Specifie the tresshold for Diferential expressed genes p-value adjusted")
    
    #Get de arguments parsed 
    args= parser.parse_args() 
    
    return args


def main():
    """
    Function to control the program execution
    """  
    #obtain the arguments parsed 
    arguments=parser() 
    
    matrix_path=arguments.matrix
    dess_path=arguments.dess
    plots=arguments.plots 
    prefix=arguments.pref 
    log_fold=arguments.logF 
    p_value=arguments.padj
    
    #verify that the usera has provided a count matrix 
    try: 
        count_matrix_df=pd.read_csv(matrix_path, sep="\t",header=0,index_col=0)
    except:
        raise ValueError("No count matrix specified") 
    
    if dess_path: 
        dess_matrix_df=pd.read_csv(dess_path, sep="\t",header=0, index_col=0)

    #obtain both matrixes for making the de analysis 
    count_trans , design_matrix = ut.pros_matrix(count_matrix_df, dess_matrix_df) 
    
    ut.write("Making the diferential expression analysis with pyDEseq2", "DE_analysis.log","../results/DE") 
    
    #make de de matrix with the diferential expression data
    DE_matrix, norm_matrix=ut.py_DESEQ2(count_trans,design_matrix) 
    
    #Plot the results 
    diff_exp_genes=ut.run_all_analyses(DE_matrix,norm_matrix,prefix=prefix,create=plots,pval_threshold=p_value,lfc_threshold=log_fold)
    
    #Finally we save the matrix with the diferentail expression data 
    DE_matrix.to_csv("../results/DE/DE_matrix.tsv",sep='\t', index=True, header=True)
    norm_matrix.to_csv("../results/DE/norm_matrix.tsv",sep='\t', index=True, header=True)
    diff_exp_genes.to_csv("../results/DE/significant_DE_genes.tsv",sep='\t', index=True, header=True) 
    
    #all the data search is used to recover the genomic postions of the DE genes 
    try:
        subprocess.run(['bash', 'recover_DE_seqs.sh', "../results/DE/significant_DE_genes.tsv", "/export/space3/users/vjimenez/Genomes/BdivergensROUEN_87/braker.gff3"], check=True) #Check
    except Exception as e: 
        ut.write(f"The position recovery was not made there was an error {e}", "DE_analysis.log","../results/DE")
if __name__ == "__main__":
	main()