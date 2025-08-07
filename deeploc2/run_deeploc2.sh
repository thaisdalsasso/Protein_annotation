#!/bin/bash

protein_dir="/Users/dalsasso/Desktop/Posdoc/CAU/data/References/secretomes/"
deeploc2_output_dir="/Users/dalsasso/Desktop/Posdoc/CAU/annotations/deeploc2"

for protein_file in "$protein_dir"*secretome.fa; do
   
    base_name=$(basename "$protein_file" | cut -d. -f1)

    /Users/dalsasso/Library/Python/3.9/bin/deeploc2 -f "$protein_file" -m Fast -o "${deeploc2_output_dir}${base_name}.secretome.deeploc2"

done

