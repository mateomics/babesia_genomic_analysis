#!/usr/bin/env bash 
set -e # Only to ensure scrpt executions
set -u # To avoid undefined variables usage
set -o pipefail # To avoid failed runs


#Make the preprocesing of the raw fastq files 
#Arguments: 
#      $1: directory path with the SRR directories files
#      $2: optional argument, specifies the output dir 
#The script generate one directory with all the clean SRR files 

#The uniq argument specified is the path 
if [[ $# != 1 ]];then   
    echo "The only need 1 argument, there were less or more, please try again"
    exit 1 

path=$1

#make the output directory if it does not exit yet 
if [[  $# > 2 ]]; then 
    output_clean=$2 
else 
    output_clean="$path"/../../results/clean
fi 

mkdir -p $output_clean

log="$output_clean"/clean.log

touch "$log"


#iterate over all the directories with SRR files in the path that is the reason of "/" at the final of the regular expresion 
for SRR_dir in "$path"/*/; do 
       echo "Cleaning the samples of the SRR file ${SRR_dir##/*}" >> $log
       #make sure that the next for cicle does not make empty files with the regular pattern
       shopt -s nullglob
       #Now we have to acces ever SRR file in the directory 
       for SRR in "$SRR_dir"/*.fastq.gz; do 
              #check if the file is no the paired file from another  
              if [[ "$SRR" != *_2.fastq.gz ]];then
                    #storage the path of the file number 1 or single end file
                     file_1=$SRR
                     #make the name of the output file 
                     output_1="$(basename "${file_1/.fastq.gz/_clean.fastq.gz}")" #basename is used for only storage the file name and not all the path 
                     #we check if the file name has a paired file 
                     if [[ "$SRR" == *_1.fastq.gz  ]];then
                            #if it has it we do the same but with the _2 termination
                            file_2="${SRR/_1.fastq.gz/_2.fastq.gz}" 
                            output_2="$(basename "${file_2/_2.fastq.gz/_2_clean.fastq.gz}")"
                            #make the cut adapt for paired files 
                            echo "Cleaning the paired end files ${SRR##/*} outputs: $output_1 $output_2" >> $log
                            nohup cutadapt -m 20 -q 25 -o "$output_clean"/"$output_1" -p "$output_clean"/"$output_2" $file_1 $file_2& #explorar la idea de usar qsub con un script externo como vero 
                     else 
                            echo "Cleaning the single end file ${SRR##/*} output: $output_1" >> $log
                            #also we run the cut adapt if the file were single end  
                            nohup cutadapt -m 20 -q 25 -o "$output_clean"/"$output_1" $file_1&
                     fi
              fi
       done
done 
