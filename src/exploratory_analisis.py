"""
The objetive of this python script is make the exploration of the information store in the count matrix, this in order to review the count distribution of the diferents samples
and classified the data into conditions instead of SRR ids 

"""
from sklearn.decomposition import PCA 
import seaborn as sns 
import matplotlib.pyplot as plt  
import numpy as np 
import pandas as pd 
import scanpy as sc
import argparse, os


def write(text:str, write:str="exploratory.log")->None: 
    """
    Function that write text to a specified file, log file is the default
    
    Arguments: 
        -text: the information to be written
        -write: the file that would be written 
    
    """ 
    #make sure the path already exist
    direc= os.path.join("output")
    if not os.path.exists(direc):
        os.mkdir(direc) 
    
    #get the path of the file 
    summaryfile=os.path.join(direc,write)
    
    #finaly we write the content we want 
    with open (summaryfile, "a") as sum: 
        sum.write(text + "\n")

def save_plt(plt_ob, name:str)->None:
    """
    Function for saving the plots in the script in a results plt directory
    Arguments
        -plot_ob(plt.figure): a matplot object with the figure that would be saved
        -name(string): the name for the plot file 
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
        write(f"Saving {name} plot in {plotfile}")
    except: 
        write(f"Saving {name} plot in {plotfile}")

def box_plot(data)->None: 
    """
    Function that makes the boxplot figure 
    Arguments: 
        -data(data_frame):the count matrix melted or in "large format"
    """
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
    
    
    
def density_plot(data)->None: 
    """
    Function that makes the density plot of the distribution of counts in both conditions
    Arguments: 
        -data(data_frame):the count matrix melted or in "large format"
    """ 
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
    
    
def PCA_plot(matrix:pd.DataFrame)->None:
    """
    Function that make a PCA plot 
    Arguments: 
        -matrix(data_frame): the normalized count matrix
    """
    
    matrix_t=matrix.T
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
    
    
def normalization(count_matrix:pd.DataFrame, info_matrix:pd.DataFrame) -> pd.DataFrame:
    """
    Function that make the normaliztion of the data for viewing how the samples in order with the conditions agroups
    Arguments: 
        -count_matrix(data_frame): the clean count matrix to be normalizated
        -info_matrix(data_frame): if it is not none it has the information to make a batch effect processing 
    Returns: 
        -counts_norm: the matrix normalized
    """ 
    
    #we make a fast normalization 
    if info_matrix: 
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

def gen_matrix(matrixes:list) -> pd.DataFrame: 
    """
    Function that help us to parsed the matrixes provided by the user
    """
    #verify id the matix with the srr info was specified 
    if len(matrixes) > 1: 
        raw_matrix_df=matrixes[1]
        raw_info=matrixes[0]
        batch=True
    else:  
        raw_matrix_df=matrixes[0]
        write("No where bacth file specified") 
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
    
    if batch :
        count_matrix=raw_matrix_df[columns[1][2:]]
        info_matrix=pd.DataFrame(raw_info["LibraryName"]) 
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
    else: 
        matrix_list=[raw_matrix_df]
        
    #we parse the matrixes: 
    write(f"Cleaning the matrixes with the condition names")
    count_matrix , info_matrix= gen_matrix(matrix_list) 
    
    #now is time for normalizating the data
    write("Making the normalization of the matrixes")
    norm_matrix=normalization(count_matrix, info_matrix) 
    
    if "box" in plots or "den" in plots: 
        #now we made the large format of the data this for boxplot an density graphs 
        df_logmelt= norm_matrix.melt(var_name="sample", value_name="Value") 

        #add some column for identifing the to types of replicates in the rna-seq data
        df_logmelt["Cond"]= df_logmelt["sample"].str.replace("_.$","", regex=True) 
    
    if "PCA" in plots: 
        write("Generating PCA plot")
        PCA_plot(norm_matrix)
        
    if "box" in plots:
        write("Generating box plot")
        box_plot(df_logmelt)
    if  "den" in plots:
        write("Generating density plot")
        density_plot(df_logmelt) 
        
        
    #finally we save the raw matrix that is now clean and with the conditions names: 
    output_path="../results/tables/count_matrix_clean.tsv"
    write(f"Saving the matrix in {output_path}")
    count_matrix.to_csv(output_path,sep='\t', index=True, header=True)

    
if __name__ == "__main__":
	main()