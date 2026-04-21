#!/usr/bin/bash


fasta_dir="/proteomes/Zv4anno"   
out_dir="/annotations/deeptmhmm"


tool_dir="$HOME/DeepTMHMM-Academic-License-v1.0"


for fasta_file in "$fasta_dir"/*_no_stop_codon.fa; do

    base_name=$(basename "$fasta_file" | cut -d_ -f1)

    echo "Analyzing $base_name..."

    cd "$tool_dir"

    # Run DeepTMHMM within a conda env
    conda run -n deeptmhmm python3.8 predict.py --fasta "$fasta_file" --output-dir "${out_dir}/${base_name}"
done

echo "Done!"

#conda deactivate