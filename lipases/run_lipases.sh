#!/bin/bash

for i in ../fastas4anno/*prot_renamedIDs.fa
do
    base_filename=$(basename "$i")
    hmmscan --domtblout "$base_filename.hmmscan.lipases.dom.tbl" /workspace/thais/lipases_db/dbLED.txt "$i" > "./$base_filename.hmmscan.lipases.dom.out"
    
    python3 hmmscan-parser.py "./$base_filename.hmmscan.lipases.dom.tbl" > "./$base_filename.hmmscan.lipases.parsed"

done

