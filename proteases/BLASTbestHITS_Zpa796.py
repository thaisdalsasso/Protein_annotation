#!/usr/bin/python3

import numpy as np
from Bio import SeqIO

blastinput = "Zpa796.merops.blastp" # BLAST output format: -outfmt "6 std qcovs"

besthit = {}
with open(blastinput,"r") as set1:
    for i in set1:
        i = i.rstrip()
        i = i.split("\t")
        query = i[0]
        subje = i[1]
        ident = float(i[2])
        align = float(i[3])
        startquery = int(i[6])
        endquery = int(i[7])
        evalu = float(i[10])
        bitscore = float(i[11])
        cov = float(i[12])

        if query != subje and ident >= 40.0 and cov >= 80.0:

            if query not in besthit:
                besthit[query] = i
            else:
                if bitscore >= float(besthit[query][11]):
                        besthit[query] = i

meanidentity = []
meanalignmen = []
for i in besthit.items():
    meanidentity.append(float(i[1][2]))
    meanalignmen.append(float(i[1][3]))

print ("Mean  identity:",np.mean(meanidentity))
print ("Mean align len:",np.mean(meanalignmen))

output = open(blastinput + "_id40cov80.besthit","w")
output.write("qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqcov\n")
for i in besthit.items():
    output.write("\t".join(map(str,i[1])) + "\n")

