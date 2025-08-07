#!/bin/bash


protein_dir="/home/dalsasso/data/References/secretomes/"
targetp2_output_dir="/home/dalsasso/annotations/targetp2/"

for protein_file in "$protein_dir"*secretome.fa; do
   
    base_name=$(basename "$protein_file" | cut -d. -f1)

    /home/dalsasso/softwares/targetp-2.0/bin/targetp -gff3 -format short -org non-pl -fasta "$protein_file" -prefix "${targetp_output_dir}${base_name}.secretome.targetp2"

done
