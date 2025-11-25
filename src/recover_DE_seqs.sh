#!/usr/bin/env bash 

set -e
set -u
set -o pipefail 

#Script for recover the CDS position of the DE genes with their IDs in the gff file
#Argumnets 
#   $1: path with the significants DE genes 
#   $2: path with the gff file 
#   $3: optional arguments to specify an output dir  

#if the user specified a output dir 
if [[ $# > 2 ]]; then 
    output_dir=$3
else 
    output_dir="../results/DE" 
fi 

genes_path=$1 
gff_path=$2 

out_file="$output_dir/DE_genes_pos.tsv"


#make the header of the table
echo -e "geneid\tDE_status\tchr\tstart\tend\tstrand" > $out_file

while read -r id de_s;do 
    pos=$(cut -f4,5,7,1 <(grep -w $id $gff_path))
    echo -e "$id\t$de_s\t$pos" >> $out_file;
done < <(awk 'BEGIN{FS=OFS="\t";}{if(NR>1) print $1,$NF}' $genes_path)