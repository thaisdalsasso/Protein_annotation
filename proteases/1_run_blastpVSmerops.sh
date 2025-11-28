#!/bin/bash

proteome_dir="/CAU/data/References/proteomes"
proteases_dir="/CAU/data/databases/proteases_db"

for i in $proteome_dir/*_no_stop_codon.fa
do
    base_filename=$(basename "$i" "_no_stop_codon.fa")

    # Run blastp command
    blastp -query "$i" -db $proteases_dir/merops_scan.lib -evalue 1e-4 -outfmt "6 std qcovs" -out "./$base_filename.merops.blastp"

done

