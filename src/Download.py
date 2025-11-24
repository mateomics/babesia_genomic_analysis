"""
The objective of this python script is to optimize the download of the data to the babesia_genomic_analysis project, but it still functionality for every
bioproject ID that needs to download SRR files of GSM or SRA
""" 

import argparse
from Bio import Entrez 
import pandas as pd
import utils_DW as ut

def parser():
    
    """
    Function that parse our arguments in order to have them in variable utils for our code 
    returns: 
        -args: an parser object with all the parsed arguments
    
    """ 
    
    #We initialize de object parser 
    parser= argparse.ArgumentParser(description="Parser to store the user information") 

    #arguments to be specified
    parser.add_argument("--project",
                        help="Specifie the ID of the Bioproject to use the script"
                        )
    parser.add_argument("--email",
                        default="ismadls@lcg.unam.mx",
                        help="Specifie the email those wich is going to be emplay for the search")
    parser.add_argument("--dbs",
                        default="gds,sra,genome,biosample",
                        help="Specifie separated with comas the databases to obtain information about their links to the Bioproject")
    parser.add_argument("--gds",
                        default= None, 
                        type=str,
                        help="Specifie to the program of how many GSM experiments of the GSE elink related files have to download or 'all' for download all files related")
    parser.add_argument("--sra",
                        default= None,
                        type=str,
                        help="Specifie to the program how many SRR related files download from de SRA elink or 'all' for download all files related")
    parser.add_argument("--biosmp",
                        default= None,
                        type=str, 
                        help="Specifie to the program how many SRR related files download from the biosamples or 'all' for download all files related")
    parser.add_argument("--ref", 
                        action="store_true",
                        help="Flag that specifie the download of the reference genome")
    parser.add_argument("--organism", 
                        default=None,
                        type=str,
                        help="Specifie the organism name for the search quoated by quotes")
    parser.add_argument("--per_dw",
                        action="store_true",
                        help="Use the interactive mode of the script, helpfull to download specific SRR related to the data bases given")
    parser.add_argument("--py",
                        action="store_true",
                        help="Indicates that only-python strategy to download SRR files is desired. When not used, a mixed python-bash strategy will be used.")
    
    #Get de arguments parsed 
    args= parser.parse_args() 
    
    return args 

#============================================================================================================================================================================

def main(): 
    """
    The main function serve as the template for running all the task in the downaload of SRR samples files and reference genome pipeline 
    
    """  
    #Obtain all the arguments given from the user
    arguments=parser() 
    
    Entrez.email= arguments.email
    projec_input=arguments.project  
    gds=arguments.gds
    sra=arguments.sra
    biosample=arguments.biosmp
    dbs=arguments.dbs.split(sep=",")
    ref=arguments.ref
    organism=arguments.organism
    personalized=arguments.per_dw
    py_only = arguments.py
    

    #Obtain the UID of the Bioproject only the first One 
    
    Bioproject_uid=ut.searcher("bioproject",projec_input)
    
    #Verify we are working with only one UID 
    if len(Bioproject_uid) != 1: 
        Bioproject_uid=str(Bioproject_uid[0])
        ut.summary(f"\n The {projec_input} has more tha one bioprojec UID, it has been used the first one")
    else:
        ut.summary(f"\nThe uid for the {projec_input} is {Bioproject_uid}\n")


    #Obtain the summary of the bioproject (only more relevant features)
    relevant_info=['Project_Id','Project_Acc','Project_Data_Type', 'Project_Title', 'Project_Description', 'Organism_Name']
    
    Bioproject_summary=ut.information("bioproject",Bioproject_uid, relevant_info,1) 
    #Obtain the organism name registered in the NCBI
    if not organism:
        organism=Bioproject_summary["Organism_Name"]
    
    #Show the bioproject information as a Data frame while the script is running 
    summary_df=pd.DataFrame.from_dict(Bioproject_summary, orient="index")
    print(summary_df)


    #obtain the elink information  
    Bioproject_elinks=ut.linker(dbs, Bioproject_uid,1) 


    #obatain the data bases support for download 
    dbs_support=["sra","gds","biosample"]
    dbs_interest=[db for db in dbs if db in dbs_support]
    
    #consult if the user wants to do a personalized download
    if personalized:
        uids_interest=ut.interactive(Bioproject_elinks, dbs_interest)   
        for db in uids_interest.keys():
            #change the flags 
            match db:
                case "sra": 
                    sra="perso"
                case "gds":
                    gds="perso"
                case "biosample": 
                    biosample="perso"    
    
    #make a warning if the user select more than one db support for downloading for get the SRR files 
    if sra and gds and bio_sample: 
        print("Warning you have selected more than one data base for downloading")
        ut.summary("\nThere were selected more than one data base for downloading it could be duplicated SRR files\n")
    
    
    #obtain the srr uids for download for each data base that is especified 
    for db in dbs_interest:
        #use the string for eval it as a variable for the match in order to obtain the value of the flag 
        match eval(db): 
            #if the user wants to download all the uids associated 
            case "all":
                uids=Bioproject_elinks[0][db]
                srr_ids=ut.obtain_srr(db, uids)
            #if the user had uses the personal downloading 
            case "perso": 
                uids=uids_interest[db]
                srr_ids=ut.obtain_srr(db,uids)
            case None:
                continue
            #if the user has specified a number of uids for downloading each srr asociated 
            case _: 
                uids=Bioproject_elinks[0][db][:int(eval(db))]
                srr_ids=ut.obtain_srr(db,uids)
        #download the srr in order with the db 
        if srr_ids and db=="sra":
            for srr in srr_ids:
                ut.download_srr(srr, sh_usage=not py_only) # False when '--py' is True
        elif  srr_ids and db=="gds":
            for srr in srr_ids.keys():
                ut.download_srr(srr_ids[srr], srr, sh_usage=not py_only)
        elif srr_ids and db=="biosample": 
            for bio_sample in srr_ids.keys():
                ut.summary(f"\nIt had downloaded the following SRR files related with {bio_sample}:\n")
                for srr in srr_ids[bio_sample]:
                    ut.summary(f"\nDownloaded the SRR files realted with {srr} SRA uid\n")
                    ut.download_srr(srr, sh_usage=not py_only)
    
    #make the tsv from the srr information 
    df_list=[]
    for uid in srr_ids: 
        consult_df=ut.fetcher(uid,"sra")
        df_list.append(consult_df)
    all_srrs_df=pd.concat(df_list, ignore_index=True) 
    ut.summary(ut.make_tsv(all_srrs_df), "srr_info.tsv")
           
    if ref:
        #now its time to download the referencie genome 
        if not organism:#If the organism name were no specified 
            #get the biosample uids linked with the Bioproject 
            biosample_elinks=ut.linker(["biosample"],Bioproject_uid)
            #We only select the first one 
            biosample_uid=biosample_elinks[0]["biosample"][0] 
            #Now consult the information in the database
            info_biosample=ut.information("biosample", biosample_uid)
            #Obtain the name of the organism by it first biosample 
            organism=info_biosample["Organism"]
        
        #call for the download of the reference genome 
        genome=ut.reference(organism)
        
        if not genome: 
            print(f"There was imposible download the reference genome of {organism}")
    
    print("The program has finished sucesfully")
    
    return
    
if __name__ == "__main__":
	main()