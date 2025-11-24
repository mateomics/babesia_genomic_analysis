#!/bin/bash 
set -e # Only to ensure scrpt executions
set -u # To avoid undefined variables usage
set -o pipefail # To avoid failed runs

#Make the count matrix with the coverage tables 
#Arguments: 
#   $1: path with the .count.txt files
#   $2: optional argument for specifing the output dir 
#The script generate the matrix of counts and a file with the list of the resgitered files

# Check that a path argument was provided
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <path_to_fastq_data>"
    exit 1
fi
path="$1" ## Path where is the tables *.count.txt


#if the user specified a output dir 
if [[ $# > 2 ]]; then 
    output_dir=$2
else 
    output_dir=$1
fi 

#generate a log file 
log="$output_dir/matrixes.log"
touch $log

echo -e "Intizating the generation of the matrix \nDirectory: $path" > "$log"

#Obtain the uniq IDs an names of the genes in all the coverage tables 
for file in "$path"/*.count.txt; do
    awk -F'\t' '{
        split($9, col9, ";"); # Arreglo con los campos que hay en la novena columna
        id = ""; name = ""; # Vacíos al inicio
        for (i in col9) {
	    # Quitar el prefijo si lo hay
            if (col9[i] ~ /^ID=/) id = substr(col9[i], 4);
            if (col9[i] ~ /^Name=/) name = substr(col9[i], 6);
        }
        if (id != "") print id "\t" name;
    }' "$file"
done | sort -u > "$output_dir/all_genes_with_names.txt"

echo "Making the count matrix..." >> "$log"
out_matrix="$output_dir/final_count_matrix.tsv"

#making the header of the matrix
echo -ne "GeneID\tGeneName" > "$out_matrix"
for file in "$path"/*.count.txt; do
    sample=$(basename "$file" .count.txt)
    echo -ne "\t$sample" >> "$out_matrix"
done
echo "" >> "$out_matrix" #add a \n to the matrix 

#we search the count of all the uniqs genes with this while bucle
while IFS=$'\t' read -r gen_id name; do 
    #fisrt we print in the file the ID an name of the gene if it has both  
    echo -ne "$gen_id\t" >> "$out_matrix"

    #if does not have name we only print "-" in its name
     if [[ -z "$name" ]]; then
        echo -ne "--" >> "$out_matrix"
    else
        echo -ne "$name" >> "$out_matrix"
    fi
    #seach the id in all the coverage tables
    for file in "$path"/*.count.txt; do
        #obtain only the id of the gene
         id=$gen_id
        #search in the table only the id an obtain the column with the count of each gene
        count_gene=$(grep -w "ID=$id" $file | cut -f10 ) 
        #if the gene was coverated in the table we print on the matrix
        if [[ -n $count_gene ]]; then 
            echo -ne "\t${count_gene}" >> "$out_matrix"
        #if it was not we asume the count was 0 
        else 
            echo -ne "\t0" >> "$out_matrix"
        fi 
    done 
    #print \n in the matrix for the next gene 
    echo "" >> "$out_matrix"
done < "$path/all_genes_with_names.txt"


echo "Matrix created: $path/final_count_matrix.tsv" >> "$log"