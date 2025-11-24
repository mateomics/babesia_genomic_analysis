"""
Module with all the fuctions for making the dieferential expression analysis, this analysis beguin with the previsualization of the count matrix data
this task is make it by exploraotry analisis.py script 

========================================================
Functions used for the data previsualization 
========================================================

Variable Dictionary:
    -Global variables: None
    
    -Function write:
        -text: string to be written to file
        -file: output filename
        -out_dir: output directory path
        -direc: full directory path
        -summaryfile: complete file path
        -sum: file object for writing
    
    -Function save_plt:
        -plt_ob: matplotlib Figure object to save
        -name: filename for saving plot
        -direc: plots directory path
        -plotfile: complete plot file path
    
    -Function box_plot:
        -data: DataFrame with count matrix in melted/long format
        -box_plt: boxplot figure object
        -ax: axes object for the plot
        -e: exception object for error handling
    
    -Function density_plot:
        -data: DataFrame with count matrix in melted/long format
        -den_plt: density plot figure object
        -ax: individual axes object in subplot
        -con: current condition in loop
        -conditions: list of experimental conditions
        -e: exception object for error handling
    
    -Function PCA_plot:
        -matrix: normalized count matrix DataFrame
        -matrix_t: transposed count matrix
        -pca: PCA object from sklearn
        -principal_components: PCA transformed data
        -Condiciones: processed condition labels
        -pca_df: DataFrame with PCA results
        -PCA_plt: PCA figure object
        -ax: axes object for PCA plot
        -e: exception object for error handling
    
    -Function normalization:
        -count_matrix: clean count matrix to normalize
        -info_matrix: metadata for batch effect processing
        -corr: flag for batch correction
        -counts_cpm: counts per million normalized data
        -counts_log: log2 transformed data
        -adata: AnnData object for batch correction
        -counts_batched_df: batch-corrected DataFrame
        -counts_norm: final normalized counts
    
    -Function gen_matrix:
        -matrixes: list of input matrix files
        -raw_matrix_df: raw count matrix DataFrame
        -raw_info: raw metadata DataFrame
        -batch: flag indicating batch correction availability
        -condition_dic: dictionary mapping conditions and batches
        -columns: list of column names from matrices
        -cols: column names from current matrix
        -raws_names: row names for indexing
        -info_matrix: processed metadata DataFrame
        -count_matrix: processed count matrix
        -col_names: generated column names for count matrix


""" 

from sklearn.decomposition import PCA 
import seaborn as sns 
import matplotlib.pyplot as plt  
import numpy as np 
import pandas as pd 
import scanpy as sc
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pathlib import Path
import argparse, os, matplotlib

"""
========================================================
Functions used for the data previsualization 
========================================================

"""

#====================================================================================================================================================================

def write(text:str, file:str, out_dir:str)->None: 
    """
    Function that write text to a specified file, log file is the default
    Parameters
    ----------
    -text:str
        The string to be writen 
    -file:str 
        The output file 
    -out_dir:str 
        The directory for the output file 
        
    Returns
    -------
    None
    
    Raises
    ------    
    """ 
    #make sure the path already exist
    direc= os.path.join(out_dir)
    if not os.path.exists(direc):
        os.mkdir(direc) 
    
    #get the path of the file 
    summaryfile=os.path.join(direc,file)
    
    #finaly we write the content we want 
    with open (summaryfile, "a") as sum: 
        sum.write(text + "\n")

#====================================================================================================================================================================

def save_plt(plt_ob, name:str)->None:
    """
    Function for saving the plots in the script in a results plt directory
    Parameters
    ----------
    -plt_ob:matplotlib.Figure
        Object with the plot to be saved 
        
    Returns
    -------
    None
    
    Raises
    ------
    """ 
    #make sure the path already exist
    direc= os.path.join("../results/plts")
    if not os.path.exists(direc):
        os.mkdir(direc) 
        
    plotfile=os.path.join(direc,name) 
    
    plt_ob.savefig(plotfile,      
            dpi=300,                    
            bbox_inches="tight")
    plt.tight_layout()
    try:
        plt.close(plt)
        write(f"Saving {name} plot in {plotfile}","exploratory.log","../results/plts")
    except: 
        write(f"Saving {name} plot in {plotfile}","exploratory.log","../results/plts")

#====================================================================================================================================================================

def box_plot(data:pd.DataFrame)->None: 
    """
    Function that makes the boxplot figure 
    Parameters
    ----------
    -data:pd.DataFrame
        The count matrix melted or in "large format"
        
    Returns
    -------
    None
    
    Raises
    ------
    """
    try:
        box_plt, ax = plt.subplots(figsize=(6, 4)) 
        sns.boxplot(data= data,ax=ax, x="sample", y="Value", hue="Cond", palette={"Merozoite": "skyblue", "IntraEri":"lightgreen"})
        
        ax.set_title("Boxplot of the counts ditribuction", fontsize=14, fontweight="bold", fontstyle="italic", color="darkblue")
        ax.set_facecolor("beige")
        ax.figure.set_facecolor("mintcream")
        
        ax.legend(
            bbox_to_anchor=(1.02, 0.5), 
            loc='upper left',          
            borderaxespad=0, 
            fontsize = 8,
            title_fontsize=20,
            labelcolor="g",
            facecolor="lightblue"
        ) 
        save_plt(box_plt,"box_exploratory.png") 
    except Exception as e: 
        write(f"Error creating the box plot {e}","exploratory.log","../results/plts")
    
#====================================================================================================================================================================   
    
def density_plot(data)->None: 
    """
    Function that makes the density plot of the distribution of counts in both conditions
    Parameters
    ----------
    -data:pd.DataFrame
        The count matrix melted or in "large format"
        
    Returns
    -------
    None
    
    Raises
    ------ 
    """ 
    try:
        den_plt=sns.displot( data=data, x="Value", hue="sample", col= "Cond", kind="kde", fill=True, alpha= 0.1)
        den_plt.fig.set_facecolor("mintcream")

        conditions=["Merozoite", "Intra eritrocyte"]
        # Personalizar cada subplot
        for ax,con in zip(den_plt.axes.flat, conditions):
            ax.set_facecolor("beige")
            ax.set_title(f"Density in {con}", fontsize=12, fontweight="bold", color="darkblue")

        # Título general
        den_plt.fig.suptitle("Distribution of Expression Counts", fontsize=16, fontweight="bold", fontstyle="italic", color="darkblue")

        den_plt._legend.set_title("Samples")
        den_plt._legend.get_frame().set_facecolor("lightblue")
        den_plt._legend.get_frame().set_edgecolor("black")
        den_plt._legend.set_bbox_to_anchor((1.15, 0.5))
        
        save_plt(den_plt,"density_exploratory.png")
    except Exception as e: 
        write(f"Error creating the density plot {e}","exploratory.log","../results/plts")

#====================================================================================================================================================================
    
def PCA_plot(matrix:pd.DataFrame)->None:
    """
    Function that make a PCA plot 
    
    Parameters
    ----------
    -matrix:pd.DataFrame
        The normalized count matrix
        
    Returns
    -------
    None
    
    Raises
    ------ 
    """
    
    matrix_t=matrix.T
    try:
        #PCA object
        pca= PCA(n_components=2) 

        principal_components=pca.fit_transform(matrix_t)

        #data frame with the 2 principal components 
        Condiciones=matrix.columns.str.replace("_.$","", regex=True) 

        pca_df=pd.DataFrame(data=principal_components, columns=["PC1", "PC2"], index=matrix_t.index)
        pca_df["Condiciones"]=Condiciones 
        
        PCA_plt, ax = plt.subplots(figsize=(6, 4))
        #We make a scaterplot for viewing the principal components 1 and 2 that explain the higher amount of variation 
        sns.scatterplot(data=pca_df, 
                        ax=ax,
                        x="PC1", 
                        y="PC2", 
                        hue="Condiciones", 
                        style="Condiciones",
                        markers={"IntraEri": "^", "Merozoite": "D"}
                        ) 

        ax.set_title("PCA plot of the varition in count matrix data", fontsize=14, fontweight="bold", fontstyle="italic", color="darkblue")
        ax.set_facecolor("beige")
        ax.figure.set_facecolor("mintcream")
        
        ax.legend(
            bbox_to_anchor=(1.02, 0.5), 
            loc='upper left',          
            borderaxespad=0, 
            fontsize = 8,
            title_fontsize=20,
            labelcolor="g",
            facecolor="lightblue"
        ) 
        save_plt(PCA_plt,"PCA_exploratory.png") 
    except Exception as e: 
        write(f"Error making the PCA {e}","exploratory.log","../results/plts")

#====================================================================================================================================================================     
    
def normalization(count_matrix:pd.DataFrame, info_matrix:pd.DataFrame, corr:int) -> pd.DataFrame:
    """
    Function that make the normalization of the data for viewing how the samples in order with the conditions groups
    Parameters
    ----------
    -count_matrix:pd.DataFrame 
        The clean count matrix to be normalizated
    -info_matrix:pd.DataFrame
        If it is not none it has the information to make a batch effect processing
        
    Returns
    -------
    -counts_norm:pd.DataFrame 
        The counts matrix normalized
        
    Raises
    ------
    """ 
    
    #we make a fast normalization 
    if corr: 
        #make a normaliztion with counts per million
        counts_cpm = (count_matrix / count_matrix.sum()) * 1e6
        counts_log = np.log2(counts_cpm + 1) 
        #create the object for making the batch correction
        adata = sc.AnnData(X=counts_log.T, obs=info_matrix) 
        #make the correction of the batch error
        sc.pp.combat(adata, key="batch") 
        #make the corrected data a Data frame
        counts_batched_df = pd.DataFrame(
            adata.X.T,                   
            index=count_matrix.index,          
            columns=count_matrix.columns       
        ) 
        counts_norm=counts_batched_df 
    #if we do not have matrix information we make a fast normalization
    else: 
        counts_norm=np.log2(count_matrix + 1) 
        
    return counts_norm
#====================================================================================================================================================================

def gen_matrix(matrixes:list) -> pd.DataFrame: 
    """
    Function that help us to parsed the matrixes provided by the user
    Parameters
    ----------
    -matrixes:list 
        List with the matrixes provided by the user
        
    Returns
    -------
    -count_matrix: pd.DataFrame
        data frame with the columns as the experimental conditions
    -info_matrix: pd.DataFrame 
        data frame with the meta data of the experimental conditions 
    
    Raises
    ------
    """
    #verify id the matix with the srr info was specified 
    if len(matrixes) > 1: 
        raw_matrix_df=matrixes[1]
        raw_info=matrixes[0]
        batch=True
    else:  
        raw_matrix_df=matrixes[0]
        write("No where bacth file specified","exploratory.log","../results/plts") 
        batch=False
    
    #dictionary with the conditions provided in the files 
    condition_dic={
        "ML1":["Merozoite_1",1],
        "ML2":["Merozoite_2",1], 
        "ML3":["Merozoite_3",2],
        "PI1":["IntraEri_1",1],
        "PI2":["IntraEri_2",1], 
        "PI3":["IntraEri_3",2]
    }
    columns=[]
    
    for matrix in matrixes: 
        #obtain the first column for making indexes
        cols=list(matrix.columns)
        raws_names=list(matrix[cols[0]])   
        matrix.index=raws_names
        #save de matrix col names 
        columns.append(cols)
    
    try:
        info_matrix=pd.DataFrame(raw_info["LibraryName"]) 
    except: 
        write("The metadat file has not the experimental conditions, it has been used the default conditions","exploratory.log","../results/plts")
        info_matrix=None
        batch=None
    
    if batch :
        count_matrix=raw_matrix_df[columns[1][2:]]

        info_matrix["Sample"]=[condition_dic[name][0] for name in list(info_matrix["LibraryName"])]
        info_matrix["batch"]=[condition_dic[name][1] for name in list(info_matrix["LibraryName"])]
        info_matrix["condition"]=info_matrix["Sample"].str.replace("_.$","",regex=True)
        info_matrix=info_matrix[["Sample","condition","batch"]]
        #make the clean count matrix with the conditions of the information matrix
        count_matrix.columns=[info_matrix.loc[srr,"Sample"] for srr in columns[1][2:]]
        info_matrix.index=list(info_matrix["Sample"])
        info_matrix.drop("Sample", axis=1, inplace=True)
    else: 
        count_matrix=raw_matrix_df[columns[0][2:]]
        col_names=["Merozoite_"+str(i) for i in range(1,4)] + ["IntraEri_" + str(i) for i in range(1,4)]
        count_matrix.columns=col_names 
        info_matrix=None 
    
    return count_matrix, info_matrix 

