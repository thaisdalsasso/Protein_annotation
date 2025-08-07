#!/bin/bash

for subfolder in */; do

    base_name=$(basename "$subfolder")

    gff_file="${subfolder}output.gff3"
    output_file="${base_name}_SP_signalp6.ids"

    if [ -e "$gff_file" ]; then
        awk '!/^#/ {print $1}' "$gff_file" > "$output_file"
        echo "Protein IDs with signal peptide in $base_name have been saved to $output_file"
    else
        echo "Error: $gff_file not found in $base_name"
    fi
done

