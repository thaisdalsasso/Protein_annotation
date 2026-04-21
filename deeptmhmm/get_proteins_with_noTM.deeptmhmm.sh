#!/bin/bash


for subfolder in */; do

    base_name=$(basename "$subfolder")

    gff_file="${subfolder}TMRs.gff3"
    output_file="${base_name}.noTM.deeptmhmm.ids"

    if [ -e "$gff_file" ]; then
        awk '/Number of predicted TMRs: 0/ {print $2}' "$gff_file" > "$output_file"
    else
        echo "Error: $gff_file not found in $base_name"
    fi
done

