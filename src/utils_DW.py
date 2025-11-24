"""
All the functions neccesary for the download of the SRR files an reference genome, it has to be used by the script Download.py
Variable Dictionary:
-Global variables: None
    
    -Function summary:
        -text: information to be written to file
        -write: filename for summary output
        -direc: directory path for output
        -summaryfile: full path to summary file
        -sum: file object for writing
    
    -Function make_tsv:
        -info: dictionary or dataframe with information to format as table
        -str_tsv: string storing the table format
        -columnas: list of column names
        -rows: list of row indices
        -info_row: list of row data
        -feature: individual data element in row
    
    -Function fetcher:
        -uid: unique identifier for the fetch
        -db: database for the fetching
        -mode: format to get from the consult
        -handle: Entrez efetch handle
        -runinfo: raw run information
        -consult_df: dataframe with consult information
        -root: XML root element
        -attr: XML attribute element
        -tag: XML tag name
        -value: XML tag value
    
    -Function information:
        -db: database for the search
        -id: ID for the search
        -info: list of fields to extract
        -show: flag to write to summary
        -handle: Entrez esummary handle
        -proj_summary: raw query results
        -consult_dir: processed query results
        -summary_df: dataframe of results
        -tsv: formatted table string
        -key: dictionary key in consult
    
    -Function searcher:
        -data_base: target database name
        -search_term: query term/ID
        -handle: Entrez search handle
        -search_results: raw search results
        -ID_search: list of found IDs
        -number_IDs: count of found IDs
    
    -Function linker:
        -dbs: list of target databases
        -id_uniq: source database ID
        -show: flag for summary writing
        -db_origin: source database name
        -dbs_summarys: summary text
        -dbs_uids: dictionary of linked IDs
        -data_base: current database in loop
        -handle: Entrez elink handle
        -linksBioProj: raw link results
        -IDs: list of linked IDs
        -only_ids: extracted ID list
        -link: individual link dictionary
    
    -Function prefetch_and_conversion:
        -srr_ids: SRR ID list to download
        -output_dir: download directory path
        -sh_usage: flag for bash vs python strategy
        -srr_id: current SRR ID in loop
        -srr_output_dir: SRR-specific directory
    
    -Function concatenate:
        -files: list of FastQ files to concatenate
        -out: output file path
        -sh_usage: flag for bash vs python strategy
        -file_out: output file object
        -p1, p2: subprocess objects
    
    -Function download_srr:
        -sra_id: list of SRA IDs to download
        -concatenate_sample: flag for file concatenation
        -sh_usage: flag for bash vs python strategy
        -handle: Entrez efetch handle
        -bite_file: raw file data
        -df: dataframe of run info
        -srr_ids: list of SRR identifiers
        -output_dir: download directory
        -files_fw: list of forward read files
        -files_rv: list of reverse read files
        -out_fw: concatenated forward output path
        -out_rv: concatenated reverse output path
        -path: base path for concatenated files
    
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
    
    -Function interactive:
        -dbs_elinks: list with database links information
        -dbs: list of allowed databases
        -dbs_dw: user-selected databases for download
        -uids_interest_all: dictionary of selected UIDs
        -db: current database in loop
        -db_uids: UIDs for current database
        -uids_interest: list of user-selected UIDs
        -state: user input loop control
        -operation: user menu choice
        -uid: specific ID for query
        -uid_summary: query results for specific ID
    
    -Function obtain_srr:
        -db: database of UIDs to link with SRR
        -db_uids: list of UIDs to consult
        -relevant_info: fields to extract from GDS
        -uid_dir: dictionary mapping UIDs to SRR IDs
        -uid: current UID in loop
        -consult_dir: query results for UID
        -samples_ids: sample accessions from GDS
        -sample: current sample in loop
        -record: SRA search results
        -bio_sample_elink: biosample to SRA links
        -srr_uids: SRR IDs from biosample
    
    -Function parser:
        -parser: argument parser object
        -args: parsed arguments
""" 

import Bio, os, urllib.request
from Bio import Entrez 
import pandas as pd
import subprocess
from io import StringIO
from ftplib import FTP
from urllib.parse import urlparse
import xml.etree.ElementTree as ET

#=============================================================================================================================================================

def summary(text:str, write:str="summary.txt")->None: 
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

def make_tsv(info:dict | pd.DataFrame)->str:
    """
    Function that makes a string with the format of a table to be write in a summary.txt 
    Arguments: 
        -dict (dict): dictionary with the information to make the table 
    returns: 
        -str_tsv (str): string with the dicctionary as a table  
    
    """
    str_tsv=""
    
    #check if the object info is a data frame
    if isinstance(info, pd.DataFrame):
        #obtain the columns
        columnas=list(info.columns)
        #make the header of the tsv with teh columns        
        for col in columnas: 
            str_tsv+=col + "\t" 
        str_tsv+="\n"
        #obtain all the rows in the data frame
        rows=list(info.index)
        #iterate over al the rows
        for row in rows:
            #obtain every column information of the row as a list 
            info_row=list(info.loc[row,:])
            #now we add this information to the tsv sting
            for feature in info_row:
                str_tsv+= feature + "\t" 
            str_tsv+="\n" 
        #return the string
        return str_tsv
    
    #if it was not a data frame we make the tsv from a dictionary
    for key in info.keys():
        str_tsv+=f"{str(key)}\t"
    
    str_tsv+="\n"
    for key in info.keys():
        str_tsv+=f"{str(info[key])}\t" 
        
    return str_tsv

#=============================================================================================================================================================

def fetcher(uid:str,db : str ="sra", mode:str = "text" )->pd.DataFrame | int : 
    """
    Function that make the efecth consult and help to obtain metadata of a specified uid in a database
    Arguments: 
        -uid: uniq identifier for the fectch
        -db(str): data base for the fetching
        -mode(str): more or formar to get from the consult
    Returns: 
        -consult_df(data_frame): data frame with the information of the single consult
        -0(int): state of a failed consult
        
    """
    
    #check if we are fetching a sra id
    if db == "sra":
        handle = Entrez.efetch(db= db, id=uid, rettype="runinfo", retmode="text")
        runinfo = handle.read()
        handle.close()

        #make sure we decode de information if it was given in binary format
        if isinstance(runinfo, bytes):
            runinfo = runinfo.decode("utf-8")
            
        #Now we read de csv with pandas
        df = pd.read_csv(StringIO(runinfo))
        #Obtain the relevant information 
        consult_df=df[['Run','Experiment','Platform','LibraryName','LibraryLayout','Sample','ScientificName','SampleName']]
        #finaly we return de data frame with the information 
        return consult_df 
    #if not we asume we are consulting xml information
    try:    
        #make the consult
        handle = Entrez.efetch(db=db, id=uid, retmode=mode)
        xml_text = handle.read()
        handle.close()
        #obtain the codified information and parse as an xml 
        root = ET.fromstring(xml_text)
        #initialice teh dictionary that will store the information 
        consult_df={}
        
        #obtain all the information in the xlm
        for attr in root.findall(".//SAMPLE_ATTRIBUTE"):
            tag = attr.findtext("TAG")
            value = attr.findtext("VALUE")
            #add the tag to the dictionary as it key 
            consult_df[tag]=list(value) 
        #make the dataframe from the dictionary 
        return pd.DataFrame.from_dict(consult_df)
        
    except: 
        return 0

#=============================================================================================================================================================
def information(db:str, id:str, info:list=[],show:int=0)->dict: 
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
      
    #Obtain the relevant fields in the consult
    if info:
        consult_dir={
            key:consult_dir.get(key) for key in info    
        }
    
    #if the consult is relevant to the summary we write it 
    if show: 
        tsv=make_tsv(consult_dir)
        summary(f"\n{tsv}\n")
    
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

def linker(dbs:list,id_uniq:str|list,show=0,db_origin="bioproject")->list: 
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
    dbs_summarys=""
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
            #Get the information as a string for the summary 
            dbs_summarys+=f"The {db_origin} {id_uniq} has {len(IDs)} uids linked with {data_base}\nthe fisrts ones are:{only_ids[:5]}\n"
            #Append all the uids in the diccionari of uids 
            dbs_uids[data_base]=only_ids
        else:
            dbs_summarys+=f"The bioporject {id_uniq} has not link with {data_base}\n"
    
    #write the result of the consult if it has to be in the summary 
    if show:
        
        summary(f"\nThe {db_origin} has this elinks:\n{dbs_summarys}")
    
    
    return [dbs_uids, dbs_summarys]

#====================================================================================================================================================================

def prefetch_and_conversion(srr_ids: list, output_dir: str='../data/SRR', sh_usage: bool=True) -> None:
    """
    Downloads each SRR file from a given ID list. Then, they're re-aranged as FastQ files.
    This is achieved either by using a mixed python-bash strategy or only python.

    Args:
        - srr_ids: SRR ID list (list)
        - output_dir: Desired output path (str)
        - sh_usage: Flag to specify which strategy to use (bool)
    """

    if sh_usage:
        for srr_id in srr_ids:
            # Run bash '.sh' script
            subprocess.run(['bash', 'download_single_srr.sh', output_dir, srr_id], check=True)
            # Write the fasta files download
            summary(f'\nThe SRRs files {srr_id} were download in {output_dir}/{srr_id}/\n')
        return None

    for srr_id in srr_ids:
        # Download the files as an SRR file
        subprocess.run(['prefetch', '--output-directory', output_dir, srr_id], check=True)
        
        # Transform the files in fastq 
        srr_output_dir = os.path.join(output_dir, srr_id)
    
        # This line obtain the fastq files separated in case being paired end 
        
        subprocess.run(['fastq-dump', srr_id, '-O', srr_output_dir,'--split-files', '--gzip'], check=True)
       
        # Write the fasta files download
        summary(f'\nThe SRRs files {srr_id} were download in {srr_output_dir}/\n') # Retrieve f-string

#====================================================================================================================================================================

def concatenate(files: list, out: str, sh_usage: bool=True) -> None:
    """
    Concatenate de-compressed FastQ files from a given file list and re-compresses the resulting file

    Args:
        - files: FastQ file's list (list)
        - out: Desired output path (str)
        - sh_usage: Flag to specify which strategy to use (bool)
    """
    if sh_usage:
        subprocess.run(['bash', 'concatenate_fastqs.sh', out, *files], check=True) # Send all files on list as separate arguments
        return None
    # Open the file ../data/SRR/sample*_*.fastq.gz
    with open(out, 'wb') as file_out: # write binary

        # Zcat give us the information about both files in the list we store, the standard_ouput (files merged)
        p1 = subprocess.Popen(['zcat'] + files, stdout=subprocess.PIPE) 
        # Now with gzip we transform the output of zcat in the file with all the runs of a sample
        p2 = subprocess.Popen(['gzip'], stdin=p1.stdout, stdout=file_out)
        # Close p1's channel
        p1.stdout.close()
        # Make sure all info processes have ended correctly
        p2.communicate()

#====================================================================================================================================================================

def download_srr(sra_id: list, concatenate_sample: bool | str =None, sh_usage: bool=True) -> None:
    """
    Function that make de download of SRR files asociated with a one or more srr uid
    
    Arguments: 
        - sra_id : list with the uid to be consulted (list)
        - concatenate: Flag to specify if files are from a determined GSM experiment. _1 denotes forward-read files.
                    Thus, _2 is used for reverse-read files (none/str)
    
    """
    
    # Use 'efetcth' to obtain the information of the UID in SRA
    handle = Entrez.efetch(db="sra", id=sra_id, rettype="runinfo", retmode="text") # Retrieve IDs info.
    bite_file= handle.read()
    handle.close()
    
    # Decode file, if binary
    if isinstance(bite_file, bytes):
        bite_file = bite_file.decode("utf-8")

    # Make a data frame with the string 
    df = pd.read_csv(StringIO(bite_file))
        
    # Visualize columns
    print(df[['Run','Experiment','Platform','LibraryName','LibraryLayout','Sample','ScientificName','SampleName']])
    
    # Obtain only SRR IDs related 
    srr_ids = list(df["Run"]) # Casting from pd.series to list, to enhance 'Bio' usage
    
    # SRR files download and conversion to FastQs
    os.makedirs(output_dir := '../data/SRR', exist_ok=True) # Makes directory not existing yet
    prefetch_and_conversion(srr_ids, output_dir, sh_usage)

    # Concatenate files if they belong to the same experiment
    if concatenate_sample: 
        # Each one of this lists conatin the two path of the fastq files download before 
        files_fw = [f'{output_dir}/{id}/{id}_1.fastq.gz' for id in srr_ids] # Como FASTA
        files_rv = [f'{output_dir}/{id}/{id}_2.fastq.gz' for id in srr_ids]

        # Create the output files with the name of the sample 
        out_fw = (path := f'{output_dir}/sample{concatenate_sample}') + '_1.fastq.gz'
        out_rv = path + '_2.fastq.gz'
    
        # Open the file ../data/SRR/sample<nameofsample>_1.fastq.gz
        concatenate(files_fw, out_fw, sh_usage)

        # And now the same with reversed-read data 
        concatenate(files_rv, out_rv, sh_usage)

        # Show concatenated fw-fw and rv-rv data
        print(f'Files concatenated in {output_dir}/sample{concatenate_sample}_<1/2>.fastq.gz')
    else: 
        print(f'SRRs files were no concatenated') 
    return None

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

def interactive(dbs_elinks:list, dbs:list)->dict:
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
                        summary(f"They are no srr asociated with {sample} GSM uid")
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
                    summary(f"They are no srr asociated with {uid} Biosample")
            return uid_dir