#!/usr/bin/env bash 
: '
Make the preprocing od the raw fastq files 
Arguments: 
    $1 — directory path with the SRR directories files

The script generate one directory with all the clean SRR files 
'

set -e
set -u
set -o pipefail  

#The uniq argument specified is the path 
if [[ $# != 1 ]];then   
    echo "The only need 1 argument, there were less or more, please try again"
    exit 1 

path=$1

#make the output directory if it does not exit yet 
[[ ! -d "$path"/clean ]] && mkdir "$path"/clean 

#storage the output path in a variable 
output_clean="$path"/clean

#iterate over all the directories with SRR files in the path that is the reason of "/" at the final of the regular expresion 
for SRR_dir in "$path"/*/; do 
       #Now we have to acces ever SRR file in the directory 
       for SRR in "$SRR_dir"/*.fastq.gz; do 
              #check if the file is no the paired file from another  
              if [[ "$SRR" != *_2.fastq.gz ]];then
                    #storage the path of the file number 1 or single end file
                     file_1=$SRR
                     #make the name of the output file 
                     output_1="$(basename "${file_1/.fastq.gz/_clean.fastq.gz}")" #basename is used for only storage the file name and not all the path 
                     #we check if the file name has a paired file 
                     if [[ -f "${SRR/_1.fastq.gz/_2.fastq.gz}"  ]];then
                            #if it has it we do the same but with the _2 termination
                            file_2="${SRR/_1.fastq.gz/_2.fastq.gz}" 
                            output_2="$(basename "${file_2/_2.fastq.gz/_2_clean.fastq.gz}")"
                            #make the cut adapt for paired files 
                            nohup cutadapt -m 20 -q 25 -o "$output_clean"/"$output_1" -p "$output_clean"/"$output_2" $file_1 $file_2& #explorar la idea de usar qsub con un script externo como vero 
                     else 
                            #also we run the cut adapt if the file were single end  
                            nohup cutadapt -m 20 -q 25 -o "$output_clean"/"$output_1" $file_1&
                     fi
              fi
       done
done