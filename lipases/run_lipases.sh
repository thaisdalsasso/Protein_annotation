#!/bin/bash


proteome_dir="/Users/dalsasso/Desktop/Posdoc/CAU/data/References/proteomes"
lipases_dir="/Users/dalsasso/Desktop/Posdoc/CAU/data/databases/lipases_db"



for i in $proteome_dir/*_no_stop_codon.fa
do
    base_filename=$(basename "$i" "_no_stop_codon.fa")
    hmmscan --domtblout "$base_filename.hmmscan.lipases.dom.tbl" $lipases_dir/dbLED.txt "$i" > "./$base_filename.hmmscan.lipases.dom.out"
    
    python3 hmmscan-parser.py "./$base_filename.hmmscan.lipases.dom.tbl" > "./$base_filename.hmmscan.lipases.parsed"

done
