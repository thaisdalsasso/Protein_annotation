#!/bin/bash

#Enter each subfolder in the current directory
for subfolder in */; do

    #Set the name of the directory as a base name
    base_name=$(basename "$subfolder")

    #Search for the file output.gff3 and extract protein IDs with signal peptide
    gff_file="${subfolder}output.gff3"
    output_file="${base_name}_SP_signalp6.ids"

    if [ -e "$gff_file" ]; then
        # Extract protein IDs with signal peptide and save to the output file
        awk '!/^#/ {print $1}' "$gff_file" > "$output_file"
        echo "Protein IDs with signal peptide in $base_name have been saved to $output_file"
    else
        echo "Error: $gff_file not found in $base_name"
    fi
done

