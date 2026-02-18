#!/bin/bash

for phobius_file in *.phobius.short.out; do

    base_name=$(basename "$phobius_file" | cut -d. -f1)

    #  Check conditions
    output_file="${base_name}.noTM.SP.phobius.ids"
    awk '!/^SEQENCE/ && $2 == 0 && $3 == "Y" {print $1}' "$phobius_file" > "$output_file"

done

