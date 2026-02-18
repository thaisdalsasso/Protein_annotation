#!/bin/bash

for tmhmm_file in *.tmhmm.out; do

    base_name=$(basename "$tmhmm_file" | cut -d. -f1)

    output_file="${base_name}.noTM.tmhmm.ids"
    awk '/Number of predicted TMHs:  0/ {print $2}' "$tmhmm_file" > "$output_file"

done

