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

"""
========================================================
Functions used for the diferential expresion analysis 
========================================================
""" 

def pros_matrix(raw_counts:pd.DataFrame, design:pd.DataFrame | None)->pd.DataFrame: 
    """
    Function that precces both matrixes previous at the DE analysis
    Parameters
    ----------
    -raw_counts:pd.DataFrame 
        The count matrix provided by the user  
    -design:pd.DataFrame | None 
        Could be the design matrix created previously in the pipeline
        
    Returns
    -------
    -count_transposed: pd.DataFrame
        The count matrix transposed and filtered by counts per milion 
    -dess_matrix: pd.DataFrame
        The design matrix with the count matrix metadata
        
    Raises
    ------
    """ 
    
    write(f"There were {raw_counts.shape[0]} genes in the provided count matrix","DE_analysis.log","../results/DE")
    
    #first we filter all the counts by counts per million (more than 5 counts per million) and have repressentation in at least 3 columns
    counts_per_milion=(raw_counts/raw_counts.sum())*1000000
    df_CountFilter= raw_counts[((counts_per_milion) >= 5).sum(axis=1) >= 3]  
    
    write(f"After filtering by\nMore than 5 counts per million\nRepresentation in at least 3 columns\nThere remain {df_CountFilter.shape[0]} genes","DE_analysis.log","../results/DE")
    
    #now we can transpose the count matrix 
    count_transposed=df_CountFilter.T
    
    if isinstance(design,pd.DataFrame): 
        #eliminate the batch column 
        dess_matrix= design.drop("batch",axis=1) 
    else: 
        #we use the base names of the columns in the count matrix
        Condiciones=df_CountFilter.columns.str.replace("_.$", "", regex=True)
        dess_matrix = pd.DataFrame(Condiciones, columns=["condition"], index=count_transposed.index) 
        
    return count_transposed, dess_matrix

#====================================================================================================================================================================

def py_DESEQ2(trasposed_count_matrix:pd.core.frame.DataFrame, metadata_states:pd.core.frame.DataFrame) -> pd.core.frame.DataFrame: 
    """
    Perform differential expression analysis using PyDESeq2.
    
    This function applies DESeq2 normalization and statistical testing to identify
    differentially expressed genes between conditions.
    
    Parameters
    ----------
    -trasposed_count_matrix : pd.DataFrame
        A transposed count matrix where rows represent samples and columns represent genes.
    -metadata_states : pd.DataFrame
        A metadata DataFrame containing sample information with a 'condition' column.
    -differential_exp_matrix_output : str
        File path where the differential expression results will be saved as a CSV file.
    
    Returns
    -------
    
    -stadistical_results_df: pd.DataFrame 
        containing differential expression results with statistical metrics
    
    Raises
    ------
    Exception
        If the DeseqDataSet creation fails or statistical analysis encounters errors.

    """

    # Creating the DeseqDataSet object
    try:
        deseq_data = DeseqDataSet(
            counts = trasposed_count_matrix,
            metadata = metadata_states,
            design_factors= "condition"
        )
    # Error manager, if the matrix was not found
    except Exception as error:
        write(f"[ERROR] No matrix found {trasposed_count_matrix} - {error}","DE_analysis.log","../results/DE")
        return None

    # Data normalization and dispersion processing
    deseq_data.deseq2()
    
    # After normalization, we use a class for stadistical analysis
    try:
        de_seq_stats = DeseqStats(deseq_data, contrast = ["condition", "Merozoite","IntraEri"])
        de_seq_stats.summary() 
    except Exception as error:
        write(f"[ERROR] The transformation form DeSeqData to DeseqStats has failed - {error}","DE_analysis.log","../results/DE")
        return None

    # At the end we use the atribute ".results_df" to save the results to a data frame

    stadistical_results_df = de_seq_stats.results_df

    # Data normalization
    normalized_counts = pd.DataFrame(deseq_data.layers['normed_counts'].T, index= deseq_data.var_names, columns=deseq_data.obs_names)
    
    return stadistical_results_df, normalized_counts

#====================================================================================================================================================================

def create_volcano_plot(plot_df, pval_threshold=0.05, lfc_threshold=1.0, figsize=(8, 6), output_file=None) -> None:
    """
    Makesvolcano plot from a given differential expression table
    
    Parameters
    ----------
    -res_df: pd.DataFrame
        DataFrame containing DESeq2 results
    -pval_threshold: float
        Adjusted p value treshold for significance
    -lfc_threshold: float
        log2 fold change treshold considered for differential expression
    -figsize: tuple
        Figure size (n, m)
    -output_file: str
        Ruta para guardar la figura (opcional)
        Figure's path (including its name) to be saved at (optional)
    
    Returns
    -------
    -matplotlib.figure.Figure or None

    Raises
    ------
    """
    try:    
       # Compute -log10(padj) managing extreme values to avoid -inf's
        plot_df['log10Neg'] = -np.log10(plot_df['padj'].clip(lower=1e-300)) # 1e-300 for every value under treshold
          
        # Color configuration
        colors = {"UP": "red", "DOWN": "forestgreen", "Non-DE": "darkgray"}
        
        # Figure
        fig, ax = plt.subplots(figsize=figsize)
        
        # Create scatter plot
        sns.scatterplot(data=plot_df, x='log2FoldChange', y='log10Neg', hue='Expression', palette=colors, ax=ax, s=20, alpha=0.7)
        
        # Add treshold lines
        ax.axhline(-np.log10(pval_threshold), color='black', linestyle='--', alpha=0.8)
        ax.axvline(lfc_threshold, color='black', linestyle='--', alpha=0.8)
        ax.axvline(-lfc_threshold, color='black', linestyle='--', alpha=0.8)
        
        # Title and label configuration
        ax.set_xlabel('log2(Fold Change)')
        ax.set_ylabel('-log10(p-value ajustado)')
        ax.set_title('Volcano Plot Differentialy Expressed Genes')
        ax.grid(True, alpha=0.3)
        ax.legend(title='Expresión')
        
        # Save if specified
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            write(f"Volcano plot saved in: {output_file}","DE_analysis.log","../results/DE")
        
        return fig
        
    except Exception as e:
        write(f"Error creating volcano plot: {e}","DE_analysis.log","../results/DE")    

#====================================================================================================================================================================
  
def create_heatmap(res_df, norm_counts, figsize=(12, 8), output_file=None) -> None:
    """
    Makes a clustered heatmap of differentially expressed genes.
    
    Parameters
    ----------
    -res_df: pd.DataFrame
        DataFrame with DESeq2 results
    -norm_counts: pd.DataFrame
        Normalized counts from DESeq2 dataset
    -pval_threshold: float
        p value threshold for significant genes (default 0.05)
    -lfc_threshold: float
        Absolute log2 fold change threshold (default 1.0)
    -figsize: tuple
        Figure size (default (12, 8))
    -output_file: str, optional
        Ruta para guardar el heatmap
        Path to save heatmap (optional)
    
    Returns
    -------
    -seaborn.matrix.ClusterGrid or None
        Clustermap function object or None if no significant genes
        (in any case saves them in the specified directory)

    Raises
    ------
    """
    try:

        # Filter significative genes
        significant_genes = res_df.index
        
        if len(significant_genes) == 0:
            write("Advertencia: No se encontraron genes significativos para el heatmap","DE_analysis.log","../results/DE")
            return None
        
        write(f"{len(significant_genes)} significant genes were found","DE_analysis.log","../results/DE")
        
        # Heatmap data
        common_genes = significant_genes.intersection(norm_counts.index)
        if len(common_genes) == 0:
            write("Error: No common genes between results and normalized counts","DE_analysis.log","../results/DE")
            return None
            
        heatmap_data = np.log1p(norm_counts.loc[common_genes])
        
        # Clustermap
        graph = sns.clustermap(data=heatmap_data, 
                figsize=figsize,
                cmap='viridis', 
                z_score=0, # Z-score per row
                method='average',
                metric='euclidean')
        
        plt.title(f'Heatmap: Differentialy Expressed Genes (n={len(common_genes)})')
        
        # Save if specified
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            write(f"Heatmap saved in: {output_file}","DE_analysis.log","../results/DE")
        
        return graph
        
    except Exception as e:
        write(f"Error. Heatmap could not be made: {e}","DE_analysis.log","../results/DE")

#====================================================================================================================================================================

def run_all_analyses(res_df, norm_matrix, output_dir="../results/DE", prefix="", create=3, pval_threshold=0.05, lfc_threshold=1.0) -> pd.DataFrame:
    """
    Executes the complete plot-analysis pipeline: volcano plot and heatmap.
    
    Parameters
    ----------
    -res_df: pd.DataFrame
        DESeq2 differential expression results
    -t_matrix: pd.DataFrame
        Transposed matrix generated, to make the heatmap
    -output_dir: str
        Output directory for figures (default: "results/")
    -prefix: str
        Prefix to save the plot(s)
    -create: int
        Flag to indicate which plots to generate (0: None, 1: Volcano, 2: Heatmap, 3: Both)
    
    Returns
    -------
    -diff_exp_genes_df: pd.DataFrame 
        Table with the significant DE genes 

    Raises
    ------
    """
    import os
    os.makedirs(output_dir, exist_ok=True)
    
    # Add '_' if its not already on prefix
    if prefix and not prefix.endswith('_'):
        prefix = f"{prefix}_"

    volcano_path = os.path.join(output_dir, f"{prefix}volcano_plot.png")
    heatmap_path = os.path.join(output_dir, f"{prefix}heatmap.png")

    plot_df = res_df.copy()
    #filter the genes 
    # Defining expression categories
    conditions = [
            (plot_df['padj'] < pval_threshold) & (plot_df['log2FoldChange'] > lfc_threshold), # UP
            (plot_df['padj'] < pval_threshold) & (plot_df['log2FoldChange'] < -lfc_threshold), # DOWN
        ]
    exp_type = ['UP', 'DOWN']
    plot_df['Expression'] = np.select(conditions, exp_type, default='Non-DE') # vectorized elif 
    
    diff_exp_genes_df=plot_df[plot_df['Expression'] != "Non-DE"] 
    
    match create:
        case 0: # None
            write("No plot saved","DE_analysis.log","../results/DE")
            return None
        
        case 1: # Volcano plot
            create_volcano_plot(plot_df,pval_threshold=pval_threshold, lfc_threshold=lfc_threshold, output_file=volcano_path)

        case 2: # Heatmap
            create_heatmap(diff_exp_genes_df, norm_matrix, output_file=heatmap_path)

        case _: # Both
            create_volcano_plot(plot_df, pval_threshold=pval_threshold, lfc_threshold=lfc_threshold, output_file=volcano_path)
            create_heatmap(diff_exp_genes_df, norm_matrix, output_file=heatmap_path)

    write(f"Plots completed. Results in: {output_dir}","DE_analysis.log","../results/DE")
    return diff_exp_genes_df