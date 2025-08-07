#!/bin/bash

protein_dir="/home/dalsasso/data/References/secretomes/"
effectorp3_output_dir="/home/dalsasso/annotations/effectorp3/"


for protein_file in "$protein_dir"*secretome.fa; do
   
    base_name=$(basename "$protein_file" | cut -d. -f1)

    python3 /home/dalsasso/softwares/EffectorP-3.0-main/EffectorP.py  -f -i "$protein_file" -o "${effectorp3_output_dir}${base_name}.secretome.effectorp3"

done






