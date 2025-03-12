#!/bin/bash

protein_dir="../../data/References/proteomes/"

for protein_file in "$protein_dir"*_no_stop_codon.fa;
do
    
    base_name=$(basename "$protein_file" | cut -d_ -f1)

    echo "Analyzing $base_name..."

    signalp6 -ff "$protein_file" -od "./${base_name}" -fmt txt -org eukarya -m fast


done

