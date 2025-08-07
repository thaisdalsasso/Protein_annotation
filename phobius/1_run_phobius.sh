#!/bin/bash

protein_dir="/home/dalsasso/data/References/proteomes/"
phobius_output_dir="/home/dalsasso/annotations/phobius/"


for protein_file in "$protein_dir"*no_stop_codon.fa; do
   
    base_name=$(basename "$protein_file" | cut -d_ -f1)

    perl /home/dalsasso/softwares/phobius/phobius.pl -short "$protein_file" > "${phobius_output_dir}${base_name}.phobius.short.out"

done


