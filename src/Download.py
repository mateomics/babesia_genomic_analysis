""""
The objetive of this python script is to optimize de download of the data to the babesia_genomic_analysis proyec, but it still functionality for ever
bioproject ID that needs to download SRR files of GSM or sra

Variable Dictionary:
    -Global variables: None
    
    -Function summary:
        -text: information to be written to file
        -write: filename for summary output
        -direc: directory path for output
        -summaryfile: full path to summary file
        -sum: file object for writing
    
    -Function information:
        -db: database name for search
        -id: ID for database query
        -info: list of fields to extract
        -show: flag to write to summary
        -handle: Entrez query handle
        -proj_summary: raw query results
        -consult_dir: processed query results
        -summary_df: dataframe of results
    
    -Function searcher:
        -data_base: target database name
        -search_term: query term/ID
        -handle: Entrez search handle
        -search_results: raw search results
        -ID_search: list of found IDs
        -number_IDs: count of found IDs
    
    -Function linker:
        -dbs: list of target databases
        -ID_uniq: source database ID
        -show: flag for summary writing
        -dbs_summarys: summary text
        -dbs_uids: dictionary of linked IDs
        -data_base: current database in loop
        -handle: Entrez elink handle
        -linksBioProj: raw link results
        -IDs: list of linked IDs
        -only_ids: extracted ID list
    
    -Function download_srr:
        -sra_id: SRA ID(s) to download
        -concatenate: flag for file concatenation
        -handle: Entrez efetch handle
        -bite_file: raw file data
        -df: dataframe of run info
        -srr_ids: list of SRR identifiers
        -output_dir: download directory
        -id: current SRR ID in loop
        -srr_output_dir: SRR-specific directory
        -files_fw: list of forward read files
        -files_rv: list of reverse read files
        -out_fw: concatenated forward output
        -out_rv: concatenated reverse output
        -file_out: output file object
        -p1, p2: subprocess objects
    
    -Function reference:
        -organism: organism name for reference
        -genome_uids: assembly database IDs
        -assembly_data: filtered assembly info
        -relevant_info: fields to extract
        -uid: current assembly ID in loop
        -uid_summary: assembly information
        -ftp_link: FTP path to genome
        -parsed: parsed FTP URL
        -ftp_server: FTP hostname
        -ftp_path: FTP directory path
        -ftp: FTP connection object
        -files: list of FTP files
        -output_dir: genome download directory
        -fasta_files: genomic FASTA files
        -gff_files: annotation GFF files
        -fasta_url: full FASTA URL
        -fasta_out: local FASTA path
        -gff_url: full GFF URL
        -gff_out: local GFF path
    
    -Function parser:
        -parser: argument parser object
        -args: parsed arguments
    
    -Function main:
        -arguments: parsed command line args
        -projec_input: bioproject ID input
        -gds: GDS download flag
        -sra: SRA download flag
        -dbs: databases to query
        -Bioproject_uid: bioproject unique ID
        -relevant_info: summary fields
        -Bioproject_summary: project information
        -organism: organism name
        -summary_df: project summary dataframe
        -Bioproject_elinks: database links
        -db_interest: current database
        -Bioprjoect_elink_GSE: GSE links
        -state: user input loop control
        -uids_interest: selected IDs
        -operation: user menu choice
        -uid: specific ID for query
        -consult_dir: query results
        -samples_ids: sample accessions
        -sample: current sample in loop
        -record: SRA search results
        -Bioprjoect_elink_sra: SRA links
        -genome: reference download success

""" 



import Bio, argparse, sys, os, urllib.request
from Bio import Entrez 
import pandas as pd
import subprocess
from io import StringIO
from ftplib import FTP
from urllib.parse import urlparse

#=============================================================================================================================================================

def summary(text:str, write="summary.txt"): 
    """
    Function to write a file with a summary of the search 
    
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
        sum.write(text)


#=============================================================================================================================================================
def information(db:str, id:str, info=[],show=0)->dict: 
    """
    Function that obtain the summary of a consult in a NCBI db according with the info to be specified
    
    Arguments: 
        -db : The data base for the search (str)
        -id : the id for the search (str)
        -info : a list with the info to be displayed, it can be empty for untrimmed consults (list) 
        -show: flag to know if the consult would be written (int)
    Returns: 
        summary_df: a data frame thah sumaraize de consult (data frame)
    
    """
    #Make the consult of the information in the database 
    try:
        handle = Entrez.esummary(db=db, id=id)
    except:
        summary(f"\nIt had been imposible to obtain the {id} summary in the {db} database\n")
        return 0
    proj_summary = Entrez.read(handle) 
    
    #Get only the dictionary with the information 
    try:
        #this for the consults of bioproject and assembly
        consult_dir=proj_summary["DocumentSummarySet"]["DocumentSummary"][0] 
    except: 
        #this for other consults 
        consult_dir=proj_summary[0]
    
    #if the consult is relevant to the summary we write it 
    if show: 
        summary(f"\n{consult_dir}\n")

    
    #Obtain the relevant fields in the consult
    if info:
        consult_dir={
            key:consult_dir.get(key) for key in info    
        }
    
    #return de diccionary with the information of the consult
    return consult_dir

#===========================================================================================================================================================

def searcher(data_base:str, search_term:str)->list: 
    """
    This functios serve as a bridge to link a raw ID to a data base specific ID 
    
    Arguments: 
        -data_base: The name of the data base to be used in the search (str)
        -search_term: the raw id to be search in the data base (str)
    Returns: 
        -ID_uniq: the ID that serves to identify that consult (list)
    """
    
    #Consult the information of the raw ID
    try:    
        handle = Entrez.esearch(db=data_base, term=search_term)
    except: 
        raise ValueError(f"The {search_term} is not yet asociated with a {data_base} UID")
    
    #Obtain the search as a python object, this to be processed 
    search_results = Entrez.read(handle)
    handle.close()  

    #Obtain the list with the ID associated 
    ID_search=search_results["IdList"]

    #Verify the ID is uniq or specify we continue with the first ID in the list 
    number_IDs=len(ID_search)
    
    #we verify it has relation with 
    if number_IDs != 1: 
        if number_IDs < 1:
            summary(F"\nThe {search_term} is no current asociated with a {data_base} uid\n")
            return 0 
        else:
            return ID_search
    else: 
        return ID_search 
    
#==============================================================================================================================================================

def linker(dbs:list,id_uniq:str,show=0,db_origin="bioproject")->list: 
    """
    Function that obtain the links with other data bases that are especified by the users 

    Arguments: 
        -dbs: list with the data bases for the consult (list)
        -ID_uniq: the id of the bioporject for search (list)
        -show: flag to know if the consult would be written (int)
        
    Returns:
        -dbs_uids: diccionary with the uids of the data bases linked with exit (dir)
    
    """
    #store the summarize of the consult (numer of uids and database)
    dbs_summarys=[]
    #store the overall uids consult as a diccionary 
    dbs_uids={}
    for data_base in dbs: 
        #make the consult 
        try:
            handle = Entrez.elink(dbfrom=db_origin, db=data_base, id=id_uniq)
        except:
            summary(f"It had been imposible make the elink between {id_uniq} bioproject and {data_base}")
        linksBioProj = Entrez.read(handle)
        handle.close()
                
        #Verifying it has link with the data base 
        if linksBioProj[0]["LinkSetDb"]:
            #Obtain de dictionary with the uids in the current data base 
            IDs=linksBioProj[0]["LinkSetDb"][0]["Link"]
            #Create a list with only the uids
            only_ids = [link["Id"] for link in IDs]
            #Get the information as a string for teh summary 
            dbs_summarys.append(f"The {db_origin} {id_uniq} has {len(IDs)} uids linked with {data_base}\nthe fisrts ones are:{only_ids[:5]}\n ")
            #Append all the uids in the diccionari of uids 
            dbs_uids[data_base]=only_ids
        else:
            dbs_summarys.append(f"The bioporject {id_uniq} has not link with {data_base}\n") 
    #write the result of the consult if it has to be in the summary 
    if show:
        summary(f"\nThe {db_origin} has this elinks:\n{dbs_summarys}")
    
    
    return [dbs_uids, dbs_summarys]
#====================================================================================================================================================================

def download_srr(sra_id:list,concatenate=None): 
    """
    Function that make de download of SRR files asociated with a one or more srr uid
    
    Arguments: 
        -sra_id : list with the uid to be consulted (list)
        -concatenate : a flag to specifie if the files are from one GSM experiment (none/str)
    
    """
    
    # Use ecfecth to obtain the information of the uid in sra
    
    handle = Entrez.efetch(db="sra", id=sra_id, rettype="runinfo", retmode="text")
    bite_file= handle.read()
    handle.close()
    
    #Its common this files are binary strings so we decodifie it 
    if isinstance(bite_file, bytes):
        bite_file = bite_file.decode("utf-8")

    #Make a data frame with the string 
    df = pd.read_csv(StringIO(bite_file))
    
        
    print(df[['Run','Experiment','Platform','LibraryName','LibraryLayout','Sample','ScientificName','SampleName']])
    
    #obtain only the SRR ids related 
    srr_ids=list(df["Run"])
    
    
    #The directory for the output 
    output_dir = "../data/SRR"

    for srr_id in srr_ids:
        #Download the files as a srr file
        subprocess.run(["prefetch", "--output-directory", output_dir, srr_id], check=True)
        
        #Transform the files in fastq 
        srr_output_dir = os.path.join(output_dir, srr_id)
    
        #This line obtain the fastq files separated in case being paired end 
        
        subprocess.run(["fastq-dump", srr_id, "-O", srr_output_dir,"--split-files", "--gzip"], check=True)
       
        
        #write the fasta files download
        summary(f"\nThe SRRs files {srr_id} were download in {srr_output_dir}/\n")

    #if the files belong to a single experiment, the we have to concatenate them 
    if concatenate: 
        #Use _1 to the files fw and _2 to the reversed files 
        #Each one of this list conatin the two path of the fastq files download before 
        files_fw = [output_dir + "/" + f + "/" + f + "_1.fastq.gz" for f in srr_ids]
        files_rv = [output_dir + "/" + f + "/" + f + "_2.fastq.gz" for f in srr_ids]

        # Create the output files with the name of the sample 
        out_fw = output_dir + f"/sample{concatenate}_1.fastq.gz"
        out_rv = output_dir + f"/sample{concatenate}_2.fastq.gz"
    
        # Open the file ../data/SRR/sample<nameofsample>_1.fastq.gz
        with open(out_fw, "wb") as file_out:
            #with zcat we obtain the information of both files in the list we store, the standar_ouput (files merged)
            p1 = subprocess.Popen(["zcat"] + files_fw, stdout=subprocess.PIPE) 
            #now with gzip we transform the output of zcat in the file with all the runs of a sample
            p2 = subprocess.Popen(["gzip"], stdin=p1.stdout, stdout=file_out)
            #Close the p1 channel
            p1.stdout.close()
            #Make sure all info being correct 
            p2.communicate()

        #We do the same with the reversed data 
        with open(out_rv, "wb") as fout:
            p1 = subprocess.Popen(["zcat"] + files_rv, stdout=subprocess.PIPE)
            p2 = subprocess.Popen(["gzip"], stdin=p1.stdout, stdout=fout)
            p1.stdout.close()
            p2.communicate()
        #now we have concatenate fw data with fw data and rv with rv in two files 
        print(f"Files concatenated in {output_dir}/sample{concatenate}_<1/2>.fastq.gz")
    else: 
        print(f"SRRs files were no concatenated") 
    return
        
        
#====================================================================================================================================================================
   
def reference(organism:str): 
    """
    Function that make the seach of the reference genome and the download if it exits 
    
    Arguments: 
        -organism: the name of the organism of the reference genome (str)
    Retunrs: 
        -1: the correct download of the genome 
        -0: there were an error in the download
    """
    
    #obatin the uids related with the genome of the bioproject organism
    genome_uids=searcher("assembly", organism )

    #list to store the refseq genomes
    assembly_data = []
    #information to the summary consult
    relevant_info=["Organism","RefSeq_category","FtpPath_GenBank"]
    #for al the uids we consult their information
    for uid in genome_uids: 
        #obtain the dictionary with the information 
        uid_summary=information("assembly", uid, relevant_info) 
        #filter only thoes which are refseq
        if uid_summary["RefSeq_category"] != "na": 
            assembly_data.append(uid_summary)


    if not assembly_data: 
        print(f"{organism} does not have a reference genome yet")
        return 0
    else: 
        #we use the fisrt refseq genome in the case it were more than one 
        ftp_link=assembly_data[0]["FtpPath_GenBank"] 
        summary(f"\nThe references genomes were:{assembly_data}\nIt has been used the first one\n")


    
    #We use urllib to parsed the link
    parsed = urlparse(ftp_link)

    #obatain the information of the parsed link  
    ftp_server = parsed.hostname        
    ftp_path = parsed.path            

    summary(f"\nThe FTP_server: {ftp_server}, with the path: {ftp_path}\n")

    # conect to the FTP server in an anonimous form 
    ftp = FTP(ftp_server)
    ftp.login() 
    ftp.cwd(ftp_path)

    #Obtain all the files in the FTP link 
    files = ftp.nlst()  
    #close the ftp link 
    ftp.quit() 
    summary(f"\nThe files stored in the ftp are:\n {files}")

    #Now we filter the files GFF an fasta (fna) does relevant for us 
    fasta_files = [f for f in files if f.endswith("_genomic.fna.gz")]
    gff_files = [f for f in files if f.endswith("_genomic.gff.gz")]

    #Use the same directory of the downloaded SRR
    output_dir = "../data/genome"

    #if the files already exist with continue to the download via urllib.request.urlretrieve
    if fasta_files or gff_files:
        #generate a directory for the genome files 
        os.makedirs(output_dir, exist_ok=True)
        #download the fasta file if it exists 
        if fasta_files:
            #Generate the link related with the fasta file
            fasta_url = ftp_link + "/" + fasta_files[0]
            print("FASTA:", fasta_url)
            #Create the name of the fasta file 
            fasta_out = os.path.join(output_dir, fasta_files[0])
            #Download fasta file 
            urllib.request.urlretrieve(fasta_url, fasta_out)
            print(f"Genome were download with exit in:{fasta_out}")
        else: 
            print(f"There were no fasta files for {organism}")
        #Use the same startegy for the GFF file    
        if gff_files:
            gff_url = ftp_link + "/" + gff_files[0]
            print("GFF:", gff_url)
            gff_out = os.path.join(output_dir,gff_files[0])
            #Download GFF file 
            urllib.request.urlretrieve(gff_url, gff_out)
            print(f"Anotation file GFF were stored in {gff_out}:")
        else:
            print(f"There were no GFF files for {organism}")  
            
    return 1
#===============================================================================================================================================================

def interactive(dbs_elinks:list, dbs:list):
    """
    Function that allows the selection of specifcs uids of data bases related with the bioproject it stores that information 
    Arguments: 
        -dbs_elinks: list with both the dictionary with the data bases and their uids associated (list)
        -dbs: list of data bases allowed for the personal download (list)
    returns: 
        -uids_interest_all: dictionary with the databases and their uids selected by the users (dict)
    """
    #Ask the users what data bases is interesed in 
    dbs_dw=input(f"You have specified this data bases {str(dbs)},\nplease insert dbs interesed in download srr files separated by a coma:").split(sep=",")
     
    uids_interest_all={}
    
    for db in dbs_dw:
        #get the uids related with the data base in the iteration
        db_uids=dbs_elinks[0][db]
        if db_uids: 
            uids_interest=[]
            state=1
            print(dbs_elinks[1])
            while state: 
                operation=int(input(f"Your consulting the {db},Insert 1 to consult especific information of one uid, 2 for download specifics SRR of a GSM uid, 3 download all, 4 for exit: "))
                match operation:
                    case 1: 
                        #Give the user the options of the GSM
                        print(f"The uids {db} associated with bioproject are: {db_uids}")
                        #request if the user is interest in one uid 
                        uid=input("Insert the uid of the GSM interesde in: ")
                        #consult those information 
                        uid_summary=information(db, uid)
                        #print to the user the information 
                        print(pd.DataFrame.from_dict(uid_summary))
                        #ask if the user want to download the id consulted
                        if int(input(f"Do you want to download the SRR related with {uid}, type 1 for yes, 0 for no: ")):
                            uids_interest.append(uid)  
                    case 2: 
                        #Give the user teh options of the GSM
                        print(f"The uids GSM associated with the bioproject are: {db_uids}")
                        #Request the udis of the GSM that are interesed in dowloading theri respective SRR samples 
                        uids_interest=input(f"Insert the uids in {db} you want to download their data, splited by a coma: ").split(sep=",")
                        state=0
                    case 3: 
                        #Asign all the uids related with 
                        uids_interest=db_uids
                        state=0
                    case 4: 
                        state=0
                    case _: 
                        print("Ivalid option, try again")
            #check the list have uids for the future srr search
            if uids_interest:
                #store the uids the user is interested in acoording to the data base 
                uids_interest_all[db]= uids_interest
        else: 
            print(f"The {db} is not linked with the bioproject")

    return uids_interest_all



#===============================================================================================================================================================

def obtain_srr(db:str, db_uids:list)->list | dict: 
    """
    Function to obtain the srr uids associated with the download support data base 
    Arguments:
        -db: database of the UIDs to be linked with the SRR uids(str)
        -db_uids: list with the UIDs to consult (list)
    Returns: 
        -uid_dir: diccionary with the multiples uids of the data base (dict)
    
    """ 
    
    
    match db:
        case "gds": 
            #relevant info for the consults
            relevant_info=['Accession',"entryType","title","summary","taxon","n_samples","FTPLink","suppFile","Samples"]
            #make a dictionary that store all the gsm uid : srr realted uids
            uid_dir={}
            for uid in db_uids: 
                #Get the information of this uid
                consult_dir=information("gds",uid,relevant_info)
                #Get only the sample uid
                samples_ids=[sample['Accession'] for sample in consult_dir.get("Samples",[] )] 
                #Now make the download of ever Sample
                for sample in samples_ids:
                    #consult the sra information 
                    record=searcher("sra",sample)
                    # Make sure it sample has an SRR associated with 
                    if not record["IdList"]:
                        summary(f"They are no srr asociated with{sample}")
                    else:  
                        #add to the dictionary
                        uid_dir[uid]=record["IdList"]
            return uid_dir
        
        
        case "sra":
            #the sra uids linked are the srr ids so the return is the same 
            return db_uids
            
        case "biosample": 
            uid_dir={}
            for uid in db_uids:
                #consult the sra asociated with the biosample 
                bio_sample_elink=linker(["sra"],uid,0,"biosample")
                srr_uids=bio_sample_elink[0]["sra"]
                if srr_uids:
                    uid_dir[uid]=srr_uids
                else: 
                    summary(f"They are no srr asociated with{uid}")
            return uid_dir
    

#===============================================================================================================================================================

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
                        help="Specifie to the program of how many GSM experiments of the GSE elink related files have to download or all for download all files related")
    parser.add_argument("--sra",
                        default= None,
                        type=str,
                        help="Specifie to the program how many SRR related files download from de SRA elink or all for download all files related")
    parser.add_argument("--biosmp",
                        default= None,
                        type=str, 
                        help="Specifie to the program how many SRR related files download from the biosamples or all for download all files related")
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
    
    #Get de arguments parsed 
    args= parser.parse_args() 
    
    return args 

#============================================================================================================================================================================

def main(): 
    """
    The main function serve as the template for running all the task in the downaload of SRR samples files and reference genome pipeline 
    
    """  
    #Obatain all the arguments given from the user
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
    

    #Obtain the uid of the Bioproject only the first One 
    
    Bioproject_uid=searcher("bioproject",projec_input)
    
    #Verify we are working with only one UID 
    if len(Bioproject_uid) != 1: 
        Bioproject_uid=Bioproject_uid[0]
        summary(f"\n The {projec_input} has more tha one bioprojec UID, it has been used the first one")
    else:
        summary(f"\nThe uid for the {projec_input} is {Bioproject_uid}\n")


    #Obtain the summary of the bioproject (only more rlevant features)
    relevant_info=['Project_Id','Project_Acc','Project_Data_Type', 'Project_Title', 'Project_Description', 'Organism_Name']
    
    Bioproject_summary=information("bioproject",Bioproject_uid, relevant_info,1) 
    #Obtain the organism name registered in the NCBI
    if not organism:
        organism=Bioproject_summary["Organism_Name"]
    
    #Show the bioproject information as a Data frame
    summary_df=pd.DataFrame.from_dict(Bioproject_summary, orient="index")
    print(summary_df)


    #obtain the elink information  
    Bioproject_elinks=linker(dbs, Bioproject_uid,1) 


    #obatain the data bases support for download 
    dbs_support=["sra","gds","biosample"]
    dbs_interest=[db for db in dbs if db in dbs_support]
    
    #consult if the user wants to do a download
    if personalized:
        uids_interest=interactive(Bioproject_elinks, dbs_interest)   
        for db in uids_interest.keys():
            #change the flags 
            match db:
                case "sra": 
                    sra="perso"
                case "gds":
                    gds="perso"
                case "biosample": 
                    bio_sample="perso"    
    

    #obatain the srr uids for download for each data base that is especified 
    for db in dbs_interest:
        #use the string for eval it as a variable for the match in order to obtain the value of the flag 
        match eval(db): 
            #if the user wants to download all the uids associated 
            case "all":
                uids=Bioproject_elinks[0][db]
                srr_ids=obtain_srr(db, uids)
            #if the user had uses the personal downloading 
            case "perso": 
                uids=uids_interest[db]
                srr_ids=obtain_srr(db,uids)
            case None:
                continue
            #if the user has specified a number of uids for downloading each srr asociated 
            case _: 
                uids=Bioproject_elinks[0][db][eval(db)]
                srr_ids=obtain_srr(db,uids)
        #download the srr in order with the db 
        if srr_ids and db=="sra":
            for srr in srr_ids:
                download_srr(srr)
        elif  srr_ids and db=="gds":
            for srr in srr_ids.keys():
                download_srr(srr_ids[srr], srr)
        elif srr_ids and db=="biosample": 
            for bio_sample in srr_ids.keys():
                for srr in bio_sample:
                    download_srr(srr)
            
    if ref:
        #now its time to download the referencie genome 
        if not organism: 
            #If the organism name were no specified
            organism=input("No organism was specified at the bioproject summary plese enter the scientific name: ")
        
        #call for the download of the reference genome 
        genome=reference(organism)
        
        if not genome: 
            print(f"There was imposible download the reference genome of {organism}")
    
    print("The program has finished sucesfully")
    
    return
    
if __name__ == "__main__":
	main()