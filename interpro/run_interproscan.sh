#!/bin/bash

# v5.54-87.0

protein_dir="/home/dalsasso/data/References/proteomes/"
interpro_output_dir="/home/dalsasso/annotations/interpro/"


for protein_file in "$protein_dir"*no_stop_codon.fa; do
   
    base_name=$(basename "$protein_file" | cut -d_ -f1)

    interproscan -cpu 6 -appl SMART,SUPERFAMILY,CDD,TIGRFAM,Pfam,Coils,Gene3D \
         -i "$protein_file" -b "${interpro_output_dir}${base_name}.interpro" -f tsv
done





