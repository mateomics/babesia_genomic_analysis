#!/bin/bash
set -e # Only to ensure scrpt executions
set -u # To avoid undefined variables usage
set -o pipefail # To avoid failed runs

#Make the alignment of the cleaned SRR files
#Arguments: 
#   $1: path with the SRR files that had been cleaned before
#   $2: path with the genome of reference 
#   $3: path with the gff file with the anotation of the genes
#   $4: optional argument for specifing an output dir 

#The script generate one directory with the bam files, and another with the coverage files 
#these directories are at the same directory level of the second argument specified

if [[ $# < 2 ]]; then
    echo "There where not ptovided all the arguments"
    exit 1
fi
path="$1"
genome=$2

#if the user specified a output dir 
if [[ -n $4 ]]; then 
    output_dir=$4 
else 
    output_dir="$path"/../aligned
fi 

mkdir -p $output_dir


log="$output_dir/alignments.log"
#make the temaporary log file for documating the prgram execution
touch "$log"


#make the uniq prefix for idetifing index genome files 
index=$(basename ${genome%%.*})_index 
#make the genome index for the aligment with the bwa tool 
bwa index -p $index $genome 

#Array for controlling the proccesed files
declare -A processed

#we search in all the fastq files
for file in "$path"/*.fastq.gz; do
    #Obtain the basename with the SRR id 
    filename=$(basename "$file" .fastq.gz)

    #Obatin only the SRR id with the extetion of clean
    base_id=$(echo "$filename" | sed 's/_1_clean$//; s/_2_clean$//; s/_clean$//; s/_1$//; s/_2$//')

    #Check if the proccesed file has been proccesed 
    if [[ ${processed["$base_id"]:-} == "yes" ]]; then
        #if it has we continue with the next 
        continue
    fi

    echo "Processing the $base_id file" >> "$log"

    #Obtain the output files name 
    file_1="${path}/${base_id}_1_clean.fastq.gz"
    file_2="${path}/${base_id}_2_clean.fastq.gz"
    file_single="${path}/${base_id}_clean.fastq.gz"

    # Verify if it is a paired end file 
    if [[ -f "$file_1" && -f "$file_2" ]]; then
        #Aling with the bwa paired end sintax
        output="${output_dir}/${base_id}.sam"
        echo "Aligning paired-end $base_id" >> "$log"
        bwa mem -o "$output" "$index" "$file_1" "$file_2"

    elif [[ -f "$file_single" ]]; then
        #If it does not have a paierd file we align with the single end bwa sintax
        output="${output_dir}/${base_id}.sam"
        echo "Aligning single-end: $base_id" >> "$log"
        bwa mem -o "$output" "$index" "$file_single"

    else
        echo "ERROR:There were not such files associted with $base_id to be aligned" >> "$log"
        continue
        fi
    #Changue the status in the array
    processed["$base_id"]="yes"
    echo "Alignment completed of: $output" >> "$log"
done

echo "All alignments completed" >> "$log"
#borramos los files con los genomas indexados 
rm -rf "$index"* 


echo -e "Starting the sam to bam conversion \nDirectory: $output_dir\n" >> "$log"

# Review all the sam files that were generated 
for sam_file in "$output_dir"/*.sam; do
    base_name="${sam_file%.*}" # Obtain the name with out the final extension

    echo "Converting the $sam_file to bam sort file" >> "$log"

    # conversion to bam
    samtools view -b "$sam_file" | samtools sort -o "${base_name}.sort.bam" && rm "$sam_file" #Make sure the old sam file been removed if the conversion was succesful
    echo "file $sam_file maked bam and sorted" >> "$log"
done 

output_tables="$output_dir"/../tables 
mkdir -p $output_tables

#check if the user pass the gff file 
if [[ -n $3 ]]; then 
    echo "GFF file were obtained, making the coveraged tables..." >> $log
    gff=$3
    for bam_file in "$output_dir"/*.sort.bam; do
        base_name=$(basename ${bam_file%%.*}) # remove all the final extension
        #make the output file for the count matrix of each bam file
        count_file="${base_name}.count.txt"
        echo "Obtaining the table $count_file of the sample $base_name" >> $log
    
        echo "Processing: $bam_file" >> "$log"
        #run the coveragedBed for the bam_file with the gff specified
        coverageBed -a "$gff" -b "$bam_file" > "$output_tables"/"$count_file"
    done 
fi

echo "Conteo de files completo" >> "$log"