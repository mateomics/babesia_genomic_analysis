"""
The objetive of this python script is make the exploration of the information store in the count matrix, this in order to review the count distribution 
of the diferents samples and classified the data into conditions instead of SRR ids

"""
import pandas as pd 
import utils_DE as ut
import argparse
             
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
                        help="Specifie the path of the count matrix")
    parser.add_argument("--batch",
                        default=None,
                        help="Specifie if do you want to make the proccesing by batch and provied the path of the srr_info.tsv")
    parser.add_argument("--plots",
                        default="PCA,box,den",
                        help="Specifie the plots do you want separated by a coma plots supported= PCA,box ,den")    
    #Get de arguments parsed 
    args= parser.parse_args() 
    
    return args

#====================================================================================================================================================================

def main():
    """
    Function to control the program execution
    """ 
    #obtain the arguments parsed 
    arguments=parser()
    
    matrix_path=arguments.matrix
    matrix_info=arguments.batch
    plots=arguments.plots.split(sep=",")
    
    #first we read both matrixes
    raw_matrix_df=pd.read_csv(matrix_path, sep="\t",header=0) 
    if matrix_info:
        raw_info=pd.read_csv(matrix_info, sep="\t", header=0)
        matrix_list=[raw_info, raw_matrix_df]
        batch_corr=1
    else: 
        matrix_list=[raw_matrix_df]
        batch_corr=0
        
    #we parse the matrixes: 
    ut.write(f"Cleaning the matrixes with the condition names","exploratory.log","../results/plts")
    count_matrix , info_matrix= ut.gen_matrix(matrix_list) 
    
    #now is time for normalizating the data
    ut.write("Making the normalization of the matrixes","exploratory.log","../results/plts")
    norm_matrix=ut.normalization(count_matrix, info_matrix, batch_corr) 
    
    if "box" in plots or "den" in plots: 
        #now we made the large format of the data this for boxplot an density graphs 
        df_logmelt= norm_matrix.melt(var_name="sample", value_name="Value") 

        #add some column for identifing the to types of replicates in the rna-seq data
        df_logmelt["Cond"]= df_logmelt["sample"].str.replace("_.$","", regex=True) 
    
    if "PCA" in plots: 
        ut.write("Generating PCA plot","exploratory.log","../results/plts")
        ut.PCA_plot(norm_matrix)
        
    if "box" in plots:
        ut.write("Generating box plot","exploratory.log","../results/plts")
        ut.box_plot(df_logmelt)
    if  "den" in plots:
        ut.write("Generating density plot","exploratory.log","../results/plts")
        ut.density_plot(df_logmelt) 
        
        
    #finally we save the raw matrix that is now clean and with the conditions names: 
    output_path="../results/tables/"

    ut.write(f"Saving the matrix in {output_path}","exploratory.log","../results/plts")
    count_matrix.to_csv(output_path+"count_matrix_clean.tsv",sep='\t', index=True, header=True)
    if batch_corr:
        info_matrix.to_csv(output_path+"design_matrix.tsv",sep='\t', index=True, header=True) 
            
if __name__ == "__main__":
	main()