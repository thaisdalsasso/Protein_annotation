#!/bin/bash

#Iterate over each *.tmhmm.out file in the current directory
for tmhmm_file in *.tmhmm.out; do

    base_name=$(basename "$tmhmm_file" | cut -d. -f1)

    #Check for the conditions
    output_file="${base_name}.noTM.tmhmm.ids"
    awk '/Number of predicted TMHs:  0/ {print $2}' "$tmhmm_file" > "$output_file"

    echo "Gene names with no transmembrane domains in $base_name have been saved to $output_file"
done

