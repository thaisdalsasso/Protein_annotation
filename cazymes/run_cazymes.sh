#!/bin/bash

proteome_dir="/CAU/data/References/proteomes"
cazymes_dir="/CAU/data/databases/cazymes_db"

for i in $proteome_dir/*_no_stop_codon.fa
do
    base_filename=$(basename "$i" "_no_stop_codon.fa")
    hmmscan --domtblout "$base_filename.hmmscan.cazymes.dom.tbl" $cazymes_dir/dbCAN.txt "$i" > "./$base_filename.hmmscan.cazymes.dom.out"
    
    python3 hmmscan-parser.py "./$base_filename.hmmscan.cazymes.dom.tbl" > "./$base_filename.hmmscan.cazymes.parsed"
done

