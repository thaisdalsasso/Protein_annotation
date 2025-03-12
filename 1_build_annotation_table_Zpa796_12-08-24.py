#!/usr/bin/python3

# Use: python3.9 script.py 

import os
from Bio import SeqIO
from bs4 import BeautifulSoup
import pandas as pd
import requests
import itertools
import csv


################################################################################
# Define input and output files:


# Functional annotation files:
signalpfile = "./signalp6/Zpa796/output.gff3"
tmhmmfile = "./tmhmm2/Zpa796.tmhmm.out"
phobiusfile = "./phobius/Zpa796.phobius.short.out"
interprofile = "./interpro/Zpa796.interpro.tsv"
targetpfile = "./targetp2/Zpa796.secretome.targetp2_summary.targetp2"
deeplocfile = "./deeploc2/Zpa796/Zpa796.secretome.deeploc2.csv"
effectorfile = "./effectorp3/Zpa796.secretome.effectorp3"
cazymefile = "./cazymes/Zpa796.hmmscan.cazymes.parsed"
lipasefile = "./lipases/Zpa796.hmmscan.lipases.parsed"
proteasefile = "./proteases/Zpa796.merops.blastp_id40cov80.besthit"
orthogroupfile = "../orthofinder/Orthogroups.tsv"  ### Update index for the species in the orthogrupsDB block of the code!!!
eggnogfile = "./eggnog/Zpa796.emapper.eggnog.tsv"
amapec_esmfold_file = "./amapec/Zpa796_esmfold_amapec_prediction.csv"
amapec_af2_file = "./amapec/Zpa796_AF2bestmodels_amapec_prediction.csv"
tmalign_file = "./tmalign/AF2_vs_ESM/Zpa796_TMscores.tsv"
antismash_dir = "./antismash7/Zpa796/"  ### Update function according the format of the gene IDs
espritz_dir = "./espritz/Zpa796/"


# Structural annotation
cath_anno_file = "./foldseek/Zpa796/easy-search/Zpa796_AF2bestmodels.vs.CATH50_foldseek_qcov-qtmscore-besthits_entry-anno.tsv"
scope_anno_file = "./foldseek/Zpa796/easy-search/Zpa796_AF2bestmodels.vs.SCOPe40_foldseek_qcov-qtmscore-besthits_entry-anno.tsv"
ecod_anno_file = "./foldseek/Zpa796/easy-search/Zpa796_AF2bestmodels.vs.ECOD40_foldseek_qcov-qtmscore-besthits_entry-anno.tsv"
pdb_anno_file = "./foldseek/Zpa796/easy-search/Zpa796_AF2bestmodels.vs.PDB_foldseek_qcov-qtmscore-besthits_pdb-anno.tsv"
af2_swiss_anno_file = "./foldseek/Zpa796/easy-search/Zpa796_AF2bestmodels.vs.AFDB-SwissProt_foldseek_qcov-qtmscore-besthits_uniprot-anno.tsv"
af2_prot_anno_file = "./foldseek/Zpa796/easy-search/Zpa796_AF2bestmodels.vs.AFDB-Proteome_foldseek_qcov-qtmscore-besthits_uniprot-anno.tsv"


# Differential expression analysis
dpi04_file= "/Users/dalsasso/Desktop/Posdoc/CAU/People/Leon_Hofmann/Zpa796/Expressed_Zpa796_DPI_04-In_Vitro_DEseq2_005.csv"
dpi07_file= "/Users/dalsasso/Desktop/Posdoc/CAU/People/Leon_Hofmann/Zpa796/Expressed_Zpa796_DPI_07-In_Vitro_DEseq2_005.csv"
dpi10_file= "/Users/dalsasso/Desktop/Posdoc/CAU/People/Leon_Hofmann/Zpa796/Expressed_Zpa796_DPI_10-In_Vitro_DEseq2_005.csv"


# Structural network
network_info_file = "/Users/dalsasso/Desktop/Posdoc/CAU/network_analysis/Zpa796/tm-score_network/pairwise_TM0.5/Zpa796_tm-scores_network_MCODE_annotations.graphml_default_node.csv"



# Fasta files:
proteome_fastafile = "../data/References/proteomes/Zpa796_no_stop_codon.fa"
secretome_fastafile = "../data/References/secretomes/Zpa796.secretome.fa"
mature_secretome_fastafile="../data/References/secretomes/mature_secretomes/Zpa796.mature_secretome.fa"

# Output table:
out_anno = open("./Zpa796_annotation_12-08-24.tsv","w") 



################################################################################
# Function to parse GenBank files for antiSMASH outputs

def parse_gbk(file_path, cluster_id):
    data = []
    for record in SeqIO.parse(file_path, "genbank"):
        for feature in record.features:
            if feature.type == "region" and "product" in feature.qualifiers:
                cluster_type = feature.qualifiers['product'][0]
                for gene_feature in record.features:
                    if gene_feature.type == "CDS":
                        gene_id = gene_feature.qualifiers.get('gene', [''])[0].split('_')[-1] # Change according gene IDs
                        gene_function = gene_feature.qualifiers.get('gene_kind', ['other'])[0]
                        data.append((gene_id, cluster_type, gene_function, cluster_id))
    return data


################################################################################
# Process Biosynthesis gene cluster (BGC) from antiSMASH annotation

antismashDB = {}
def process_directory(directory):
    for filename in os.listdir(directory):
        if filename.endswith(".gbk") and "region" in filename:
            file_path = os.path.join(directory, filename)
            cluster_id = filename.replace('.gbk', '')  # Removing the file extension to get the cluster ID
            parsed_data = parse_gbk(file_path, cluster_id)
            for gene_id, cluster_type, gene_function, cluster_id in parsed_data:
                if gene_id not in antismashDB:
                    antismashDB[gene_id] = [cluster_type, gene_function, cluster_id]
process_directory(antismash_dir)

print("Number of genes in antiSMASH clusters:", len(antismashDB))


################################################################################
################################################################################
#Parse FASTA sequences used for predictions

proteomeDB = {}
with open(proteome_fastafile,"r") as set1:
    for i in SeqIO.parse(set1, "fasta"):
        seqID = i.id
        seq_length = len(i.seq)
        seq_fasta = i.seq
        proteomeDB[seqID] = [seq_length, seq_fasta]
print ("Number of protein sequences:", len(proteomeDB))


secretomeDB = {}
with open(secretome_fastafile,"r") as set1:
    for i in SeqIO.parse(set1, "fasta"):
        seqID = i.id
        fasta = i.seq
        secretomeDB[seqID] = [i.seq]
print ("Number of secreted proteins:", len(secretomeDB))


maturesecretomeDB = {}
with open(mature_secretome_fastafile,"r") as set1:
    for i in SeqIO.parse(set1, "fasta"):
        seqID = i.id
        seq_length = len(i.seq)
        seq_fasta = i.seq
        maturesecretomeDB[seqID] = [seq_length, seq_fasta]
print ("Number of mature proteins from secretome:", len(maturesecretomeDB))


################################################################################
# Parse orthogrups

orthogrupsDB = {}
with open(orthogroupfile, "r") as set1:
    next(set1)  
    for i in set1:
        i = i.rstrip().split("\t")
        geneIDs = i[4]  # Change according to the species index in the orthogroup file
        if geneIDs != "-":  
            orthogroupID = i[0]
            geneIDs = geneIDs.split(",")
            for gene in geneIDs:
                geneID = gene.strip()  # Strip any extra whitespace
                if geneID:  
                    orthogrupsDB[geneID] = [orthogroupID]
print("Total number of genes within orthogroups:", len(orthogrupsDB))


###############################################################################
# Parse Signalp6 predictions output

signalpDB = {}
with open(signalpfile,"r") as set1:
    for i in set1:
        if "#" not in i:
            i = i.rstrip().split()
            seqID = i[0]
            SPstart= i[3]
            SPend= i[4]
            secretion = i[5]
            signalpDB[seqID] = [SPstart, SPend, secretion]
print ("Sequences with SP predicted by SignalP6:",len(signalpDB))


###############################################################################
# Parse TMHMM2 predictions output

tmhmmDB = {}
with open(tmhmmfile,"r") as set1:
    for i in set1:
        if i.startswith("#"):
            if "Number of predicted TMHs:" in i:
                j = i.rstrip().split()
                seqID = j[1]
                tmdomain = j[-1]
                tmhmmDB[seqID] = [tmdomain]
print ("Sequences with TM predicted by TMHMM2:",len(tmhmmDB))


###############################################################################
# Parse Phobius predictions output

phobiusDB = {}
with open(phobiusfile,"r") as set1:
    for i in set1:
        if "SEQENCE" not in i:
            i = i.rstrip().split()
            seqID = i[0]
            tm = i[-3]
            sp = i[-2]
            phobiusDB[seqID] = [sp, tm]
print ("Sequences with SP and/or TM annotated by Phobius:",len(phobiusDB))


###############################################################################
# Parse INTERPRO predictions output

interproDBaccessions = {} 
with open(interprofile,"r") as set1:
    for i in set1:
        i = i.rstrip().split("\t")
        sequenceID = i[0]
        ipr_accession = i[11]
        ipr_domain = i[12]
        if ipr_accession != "-":
            if sequenceID not in interproDBaccessions:
                interproDBaccessions[sequenceID] = [ipr_accession]
            else:
                interproDBaccessions[sequenceID].append(ipr_accession)
print ("Sequences with domains predicted by InterProSan:",len(interproDBaccessions))


################################################################################
# Parse TargetP2 predictions output

targetp2DB = {}
with open(targetpfile,"r") as set1:
    for i in set1:
        if "#" not in i:
            i = i.rstrip().split("\t")
            geneID = i[0]
            signal = i[1]
            targetp2DB[geneID] = [signal]
print ("Sequences annotated by TargetP2:",len(targetp2DB))


################################################################################
# Parse DeepLoc2 predictions output

deeploc2DB = {}
with open(deeplocfile,"r") as set1:
    for i in set1:
        if "Protein_ID" not in i:
            i = i.rstrip().split(",")
            geneID = i[0]
            localization = i[1]
            deeploc2DB[geneID] = [localization]
print ("Sequences annotated by DeepLoc2:",len(deeploc2DB))


################################################################################
# Parse EffectorP3 predictions output

effectorDB = {} 
with open(effectorfile,"r") as set1:
    for i in set1:
        if "#" not in i:
            i = i.rstrip().split("\t")
            sequence = i[0]
            effectorDB[sequence] = i[4]
print ("Sequences annotated by EffectorP3:",len(effectorDB))


################################################################################ 
# Parse CAZymes predictions output

cazymesDB = {} 
with open(cazymefile,"r") as set1:
    for i in set1:
        i = i.rstrip().split("\t")
        sequence = i[2]
        cazyfamily = i[0].replace(".hmm", "")
        if sequence not in cazymesDB:
            cazymesDB[sequence] = [cazyfamily]
        else:
            cazymesDB[sequence].append(cazyfamily)
print ("Sequences with domains predicted as CAZymes:",len(cazymesDB))


################################################################################
# Parse Proteases predictions output

proteasesDB = {} 
with open(proteasefile,"r") as set1:
    for i in set1:
        i = i.rstrip()
        if "qseqid" not in i:
            i = i.split("\t")
            sequence = i[0]
            blasthit = i[1]
            proteasesDB[sequence] = [blasthit]
print ("Sequences with domains predicted as proteases:",len(proteasesDB))


################################################################################
# Parse Lipases predictions output

lipasesDB = {} 
with open(lipasefile,"r") as set1:
    for i in set1:
        i = i.rstrip().split("\t")
        sequence = i[2]
        lipfamily = i[0].replace(".hmm", "")
        if sequence not in lipasesDB:
            lipasesDB[sequence] = [lipfamily]
        else: 
            lipasesDB[sequence].append(lipfamily)
print ("Sequences with domains predicted as Lipases:",len(lipasesDB))


################################################################################
# Parse eggNOG predictions output

eggnogDBfeatures = {}
with open(eggnogfile,"r") as set1:
    for i in set1:
        i = i.rstrip()
        if "#" not in i:
            i = i.split("\t")
            seqID = i[0]
            kog_categ = i[6]
            if kog_categ == "-":
                kog_categ = "NaN"
            descritption = i[7]
            if descritption == "-":
                descritption = "NaN"
            go_term = i[9]
            if go_term == "-":
                go_term = "NaN"
            kegg_path = i[12]
            if kegg_path == "-":
                kegg_path == "NaN"

            eggnogDBfeatures[seqID] = [kog_categ, go_term, descritption, kegg_path]
print ("Sequences annotated by eggNOG:",len(eggnogDBfeatures))


################################################################################
# parse AMAPEC predictions output, including pLDDT values for protein structure predicted by ESMFold and AlphaFold2

amapec_esmfoldDB = {} 
with open(amapec_esmfold_file,"r") as set1:
    for i in set1:
        if "Protein ID" not in i:
            i = i.rstrip().split(",")
            seqID = i[0]
            seq_pLDDT = i[1]
            seq_AM_prob = i[2]
            seq_AM_annot = i[3]
            amapec_esmfoldDB[seqID] = [seq_pLDDT, seq_AM_prob, seq_AM_annot]
print ("Sequences annotated by AMAPEC (ESMFold-based):",len(amapec_esmfoldDB))


amapec_af2DB = {} 
with open(amapec_af2_file,"r") as set1:
    for i in set1:
        if "Protein ID" not in i:
            i = i.rstrip().split(",")
            seqID = i[0]
            seq_pLDDT = i[1]
            seq_AM_prob = i[2]
            seq_AM_annot = i[3]
            amapec_af2DB[seqID] = [seq_pLDDT, seq_AM_prob, seq_AM_annot]
print ("Sequences annotated by AMAPEC (AlphaFold2-based):",len(amapec_af2DB))


################################################################################
# parse EspritZ predictions output

espritzDB = {}
for filename in os.listdir(espritz_dir):
    if filename.endswith(".fasta.stats"):
        protein_id = filename.replace('t1.fasta.stats', '.t1')
        disorder_rate = None
        disorder_number = None
        with open(os.path.join(espritz_dir, filename), 'r') as file:
            for line in file:
                if "Total % disorder:" in line:
                    disorder_rate = line.rstrip().split()[-1]
                if "Number of disordered segments:" in line:
                    disorder_number = line.rstrip().split()[-1]
        if disorder_rate is not None and disorder_number is not None:
            espritzDB[protein_id] = [disorder_rate, disorder_number]

print("Sequences annotated by EspritZ:", len(espritzDB))


################################################################################
# parse TM-scores and RMSD values comparing AF2 and ESMFold sctructure for each protein

tmalignDB = {} 
with open(tmalign_file,"r") as set1:
    for i in set1:
        if "Protein_Name" not in i:
            i = i.rstrip().split("\t")
            seqID = i[0]
            seq_tmscore = i[1]
            seq_rmsd = i[2]
            tmalignDB[seqID] = [seq_tmscore, seq_rmsd]
print ("Sequences annotated by TMalign:",len(tmalignDB))


################################################################################
# Retrieve entries' information for CATH50 best hits from the Foldseek structural annotation (easy search)

cathDB = {}
with open(cath_anno_file, "r") as set1:
    for i in set1:
        if "query" not in i:
            i = i.rstrip().split("\t")
            queryID = i[0].replace(".pdb", "")
            cath_entryID = i[1]
            if cath_entryID.startswith("af_"):
                cath_entryID = cath_entryID.rstrip().split("_")[1]
            cath_qtmscore = i[15]
            cath_class = i[24]
            cath_architecture = i[25]
            cath_topology = i[26]
            cath_superfamily = i[27]
            cathDB[queryID] = [cath_entryID, cath_qtmscore, cath_class, cath_architecture, cath_topology, cath_superfamily]
print ("Protein structures annoated by CATH50 with Foldseek easy-serch:", len(cathDB))


################################################################################
# Retrieve entries' information for PDB best hits from the Foldseek structural annotation (easy search)

pdbDB = {}
with open(pdb_anno_file, "r") as set1:
    for i in set1:
        if "query" not in i:
            i = i.rstrip().split("\t")
            queryID = i[0].replace(".pdb", "")
            pdb_qtmscore = i[15]
            pdb_entryID = i[24]
            pdb_description = i[25]
            pdb_classification = i[26]
            pdbDB[queryID] = [pdb_entryID, pdb_qtmscore, pdb_description, pdb_classification]
print ("Protein structures annoated by PDB with Foldseek easy-serch:", len(pdbDB))


################################################################################
# Retrieve entries' information for AF2 DB (Swiss-Prot) best hits from the Foldseek structural annotation (easy search)

af2_swissDB = {}
with open(af2_swiss_anno_file, "r") as set1:
    for i in set1:
        if "query" not in i:  
            i = i.rstrip().split("\t")
            if len(i) == 27:
                queryID = i[0].replace(".pdb", "")
                af2_swiss_qtmscore = i[15]
                af2_swiss_uniprot_entryID = i[24] 
                af2_swiss_uniprot_name = i[25]
                af2_swiss_uniprot_mol_function = i[26]
            else:
                queryID = i[0].replace(".pdb", "")
                af2_swiss_qtmscore = i[15]
                af2_swiss_uniprot_entryID = i[24] 
                af2_swiss_uniprot_name = i[25]
                af2_swiss_uniprot_mol_function = "-"
            af2_swissDB[queryID] = [af2_swiss_uniprot_entryID, af2_swiss_qtmscore, af2_swiss_uniprot_name, af2_swiss_uniprot_mol_function]
print("Protein structures annotated by AF2 DB (Swiss-Prot) with Foldseek easy-search:", len(af2_swissDB))


################################################################################
# Retrieve entries' information for AF2 DB (Proteome) best hits from the Foldseek structural annotation (easy search)

af2_proteomeDB = {}
with open(af2_prot_anno_file, "r") as set1:
    for i in set1:
        if "query" not in i:  
            i = i.rstrip().split("\t")
            if len(i) == 27:
                queryID = i[0].replace(".pdb", "")
                af2_prot_qtmscore = i[15]
                af2_prot_uniprot_entryID = i[24] 
                af2_prot_uniprot_name = i[25]
                af2_prot_uniprot_mol_function = i[26]
            else:
                queryID = i[0].replace(".pdb", "")
                af2_prot_qtmscore = i[15]
                af2_prot_uniprot_entryID = i[24] 
                af2_prot_uniprot_name = i[25]
                af2_prot_uniprot_mol_function = "-"
            af2_proteomeDB[queryID] = [af2_prot_uniprot_entryID, af2_prot_qtmscore, af2_prot_uniprot_name, af2_prot_uniprot_mol_function]
print ("Protein structures annoated by AF2 DB (Proteome) with Foldseek easy-serch:", len(af2_proteomeDB))


################################################################################
# Retrieve entries' information for SCOPe40 best hits from the Foldseek structural annotation (easy search)

scopeDB = {}
with open(scope_anno_file, "r") as set1:
    for i in set1:
        if "query" not in i:
            i = i.rstrip().split("\t")
            queryID = i[0].replace(".pdb", "")
            scope_entryID = i[1].rstrip().split(".ent")[0]
            scope_qtmscore = i[15]
            scope_class = i[23]

            if len(i) == 27:  
                scope_fold = i[24]
                scope_superfamily = i[25]
                scope_family = i[26]
            else:  
                scope_fold = "-"
                scope_superfamily = "-"
                scope_family = "-"

            scopeDB[queryID] = [scope_entryID, scope_qtmscore, scope_class, scope_fold, scope_superfamily, scope_family]
print("Protein structures annotated by SCOPe40 with Foldseek easy-search:", len(scopeDB))


################################################################################
# Retrieve entries' information for ECOD40 best hits from the Foldseek structural annotation (easy search)

ecodDB = {}
with open(ecod_anno_file, "r") as set1:
    for i in set1:
        if "query" not in i:
            i = i.rstrip().split("\t")
            queryID = i[0].replace(".pdb", "")
            ecod_qtmscore = i[15]
            ecod_entryID = i[23]
            ecod_a_group = i[24]
            ecod_x_group = i[25]
            ecod_h_group = i[26]
            ecod_t_group = i[27]
            ecod_f_group = i[28]
            ecodDB[queryID] = [ecod_entryID, ecod_qtmscore, ecod_a_group, ecod_x_group, ecod_h_group, ecod_t_group,ecod_f_group]
print("Protein structures annotated by ECOD40 with Foldseek easy-search:", len(ecodDB))


################################################################################
################################################################################
# Parse genes in sctructural clsuters tables from DESeq2 pipeline

networkDB = {}
with open(network_info_file, "r") as set1:
    reader = csv.DictReader(set1)
    for row in reader:
        gene_id = row["shared name"]
        subgraph = row["Subgraph"]
        mcode_cluster = row["MCODE::Clusters (1)"] if "Cluster" in row["MCODE::Clusters (1)"] else "NaN"
        networkDB[gene_id] = [subgraph, mcode_cluster]
print("Protein in a strcutural similarity cluster:", len(networkDB))


################################################################################
################################################################################
# Parse DE tables from DESeq2 pipeline

dpi04_expressionDB = {}
with open(dpi04_file, "r") as set1:
    for line in set1:
        if "log2FoldChange" not in i:
            i = line.strip().split(",")
            gene_id = i[0].replace('"', '') + '.t1'
            log2fc = i[2]
            pval_adj = i[-1]
        dpi04_expressionDB[gene_id] = [log2fc, pval_adj]

print(len(dpi04_expressionDB), "genes processed from file", dpi04_file)


dpi07_expressionDB = {}
with open(dpi07_file, "r") as set1:
    for line in set1:
        if "log2FoldChange" not in i:
            i = line.strip().split(",")
            gene_id = i[0].replace('"', '') + '.t1' 
            log2fc = i[2]
            pval_adj = i[-1]
        dpi07_expressionDB[gene_id] = [log2fc, pval_adj]

print(len(dpi07_expressionDB), "genes processed from file", dpi07_file)


dpi10_expressionDB = {}
with open(dpi10_file, "r") as set1:
    for line in set1:
        if "log2FoldChange" not in i:
            i = line.strip().split(",")
            gene_id = i[0].replace('"', '') + '.t1' 
            log2fc = i[2]
            pval_adj = i[-1]
        dpi10_expressionDB[gene_id] = [log2fc, pval_adj]

print(len(dpi10_expressionDB), "genes processed from file", dpi10_file)



################################################################################
################################################################################
# Write output

header = ["Protein ID", "Protein length (#AA)", "Mature protein (#AA)", "Secretome", "Disorder region (percent)", "# Disorder region(s)","pLDDT (ESMFold)", "AM probability (ESMFold/AMAPEC)", "AM prediction (ESMFold/AMAPEC)", "pLDDT (AF2)", "AM probability (AF2/AMAPEC)", "AM prediction (AF2/AMAPEC)", "TM-score (AF2 vs ESMFold)", "RMSD(AF2 vs ESMFold)","SP start", "SP end", "SP score (SignalP6)", "SP (Phobius)", "#TM (Phobius)", "#TM (TMHMM)", "InterPro domain(s)", "TargetP2", "DeeplLoc2", "EffectorP3", "CAZyme family", "Protease family", "Lipase family" ,  "BS Cluster ID (Antismash)",  "BSC function (Antismash)","Gene Function in BGC (Antismash)", "Orthogroup ID", "KOG categories", "GO terms", "KEGG pathway", "eggNOG Description", "Structural subgraph", "MCODE cluster", "CATH entry", "TM-score (CATH)", "CATH Class (code)", "CATH Architecture (code)", "CATH Topology (code)", "CATH Superfamily (code)", "SCOPe entry", "TM-score (SCOPe)", "SCOPe Class", "SCOPe Fold", "SCOPe Superfamily", "SCOPe Family", "ECOD entry", "TM-score (ECOD)", "ECOD A (Architecture)", "ECOD X (Possible Homology)", "ECOD H (Homology)", "ECOD T (Topology)", "ECOD F (Family)", "PDB entry", "TM-score (PDB)", "PDB Descritption", "PDB Classification", "Swiss-Prot (AFDB) entry", "TM-score (Swiss-Prot)", "Swiss-Prot (AFDB) Name", "Swiss-Prot (AFDB) Molecular Function", "Proteome (AFDB) entry", "TM-score (Proteome AFDB)","Proteome (AFDB) Name", "Proteome (AFDB) Molecular Function", "Log2FC (04-DPI)", "Padj (04-DPI)", "Log2FC (07-DPI)", "Padj (07-DPI)", "Log2FC (10-DPI)", "Padj (10-DPI)", "Protein sequence"]

rows = []
for geneID in proteomeDB:

    length = proteomeDB[geneID][0]
    aa_seq = proteomeDB[geneID][1]
    mature_length = ''
    orthogroup = ''
    secretome = ''
    sp_start = ''
    sp_end = ''
    signalp = ''
    tmhmm = ''
    sp2 = ''
    tm2 = ''
    interpro = ''
    targetp = ''
    deeploc = ''
    effectorp = ''
    cazymes = ''
    proteases = ''
    lipases = ''
    kog = ''
    go = ''
    desctrip = ''
    pathway = ''
    pLDDT_esmfold = ''
    prob_AMP_esmfold = ''
    annot_AMP_esmfold = ''
    pLDDT_af2 = ''
    prob_AMP_af2 = ''
    annot_AMP_af2 = ''
    tmscore = ''
    rmsd = ''
    bsc_function = ''
    bsc_gene_function = ''
    bsc_id = ''
    disorder_perc = '' 
    disorder_num = ''
    subgraph = ''
    mcode = ''
    cath_id = ''
    cath_tm = ''
    cath_class =''
    cath_archit = ''
    cath_topo = ''
    cath_family = ''
    scope_id = ''
    scope_tm = ''
    scope_class = ''
    scope_fold = ''
    scope_superfam = ''
    scope_fam = ''
    ecod_id = ''
    ecod_tm = ''
    ecod_a = ''
    ecod_x = ''
    ecod_h = ''
    ecod_t =''
    ecod_f = ''
    pdb_id = ''
    pdb_tm = '' 
    pdb_descr = ''
    pdb_classif = ''
    af2_swiss_id = ''
    af2_swiss_tm = ''
    af2_swiss_name = ''
    af2_swiss_function = ''
    af2_prot_id = ''
    af2_prot_tm = ''
    af2_prot_name = ''
    af2_prot_function = ''
    dpi04_l2fc = ''
    dpi04_padj = ''
    dpi07_l2fc = ''
    dpi07_padj = ''
    dpi10_l2fc = ''
    dpi10_padj = ''


    if geneID in maturesecretomeDB:
        mature_length = maturesecretomeDB[geneID][0]
    else:
        mature_length = "NaN"

    if geneID in secretomeDB:
        secretome = "Yes"
    else:
        secretome = "No"
    
    if geneID in orthogrupsDB:
        orthogroup = orthogrupsDB[geneID]
    else:
        orthogroup = "NaN"

    if geneID in signalpDB:
        sp_start = signalpDB[geneID][0]
        sp_end = signalpDB[geneID][1]
        signalp = signalpDB[geneID][2]
    else:
        sp_start = "NaN"
        sp_end = "NaN"
        signalp = "NaN"

    if geneID in tmhmmDB:
        tmhmm = tmhmmDB[geneID]
    else:
        tmhmm = "NaN"

    if geneID in phobiusDB:
        sp2 = phobiusDB[geneID][0]
        tm2 = phobiusDB[geneID][1]
    else:
        sp2 = "NaN"
        tm2 = "NaN"

    if geneID in interproDBaccessions:
        interpro = ','.join(interproDBaccessions[geneID])
    else:
        interpro = "NaN"

    if geneID in targetp2DB:
        targetp = ','.join(targetp2DB[geneID])
    else:
        targetp = "NaN"


    if geneID in deeploc2DB:
        deeploc = deeploc2DB[geneID]
    else:
        deeploc = "NaN"

    if geneID in effectorDB:
        effectorp = effectorDB[geneID]
    else:
        effectorp = "NaN"

    if geneID in cazymesDB:
        cazymes = ','.join(cazymesDB[geneID])
    else:
        cazymes = "NaN"

    if geneID in proteasesDB:
        proteases = ','.join(proteasesDB[geneID])
    else:
        proteases = "NaN"

    if geneID in lipasesDB:
        lipases = ','.join(lipasesDB[geneID])
    else:
        lipases = "NaN"

    if geneID in eggnogDBfeatures:
        kog = eggnogDBfeatures[geneID][0].strip('"')
        go = eggnogDBfeatures[geneID][1].strip('"')
        desctrip = eggnogDBfeatures[geneID][2].strip('"')
        pathway = eggnogDBfeatures[geneID][3].strip('"')
    else:
        kog = "NaN"
        go = "NaN"
        desctrip = "NaN"
        pathway = "NaN"

    if geneID in amapec_esmfoldDB:
        pLDDT_esmfold = amapec_esmfoldDB[geneID][0].strip('"')
        prob_AMP_esmfold = amapec_esmfoldDB[geneID][1].strip('"')
        annot_AMP_esmfold = amapec_esmfoldDB[geneID][2].strip('"')
    else:
        pLDDT_esmfold = "NaN"
        prob_AMP_esmfold = "NaN"
        annot_AMP_esmfold = "NaN"

    if geneID in amapec_af2DB:
        pLDDT_af2 = amapec_af2DB[geneID][0].strip('"')
        prob_AMP_af2 = amapec_af2DB[geneID][1].strip('"')
        annot_AMP_af2 = amapec_af2DB[geneID][2].strip('"')
    else:
        pLDDT_af2= "NaN"
        prob_AMP_af2 = "NaN"
        annot_AMP_af2 = "NaN"

    if geneID in tmalignDB:
        tmscore = tmalignDB[geneID][0].strip('"')
        rmsd = tmalignDB[geneID][1].strip('"')
    else:
        tmscore = "NaN"
        rmsd = "NaN"

    if geneID in antismashDB:
        bsc_function = antismashDB[geneID][0].strip('"')
        bsc_gene_function = antismashDB[geneID][1].strip('"')
        bsc_id = antismashDB[geneID][2].strip('"')
    else:
        bsc_function = "NaN"
        bsc_gene_function = "NaN"
        bsc_id = "NaN"

    if geneID in espritzDB:
        disorder_perc = espritzDB[geneID][0].strip('"')
        disorder_num = espritzDB[geneID][1].strip('"')
    else:
        disorder_perc = "NaN"
        disorder_num = "NaN"

    if geneID in networkDB:
        subgraph = networkDB[geneID][0].strip('"')
        mcode  = networkDB[geneID][1].strip('"')
    else:
        subgraph = "NaN"
        mcode = "NaN"

    if geneID in cathDB:
        cath_id = cathDB[geneID][0].strip('"')
        cath_tm = cathDB[geneID][1].strip('"')
        cath_class = cathDB[geneID][2].strip('"')
        cath_archit = cathDB[geneID][3].strip('"')
        cath_topo = cathDB[geneID][4].strip('"')
        cath_family = cathDB[geneID][5].strip('"')
    else:
        cath_id = "NaN"
        cath_tm = "NaN"
        cath_class = "NaN"
        cath_archit = "NaN"
        cath_topo = "NaN"
        cath_family = "NaN"

    if geneID in scopeDB:
        scope_id = scopeDB[geneID][0].strip('"')
        scope_tm = scopeDB[geneID][1].strip('"')
        scope_class = scopeDB[geneID][2].strip('"')
        scope_fold = scopeDB[geneID][3].strip('"')
        scope_superfam = scopeDB[geneID][4].strip('"')
        scope_fam = scopeDB[geneID][5].strip('"')
    else:
        scope_id = "NaN"
        scope_tm = "NaN"
        scope_class = "NaN"
        scope_fold = "NaN"
        scope_superfam = "NaN"
        scope_fam = "NaN"


    if geneID in ecodDB:
        ecod_id = ecodDB[geneID][0].strip('"')
        ecod_tm = ecodDB[geneID][1].strip('"')
        ecod_a = ecodDB[geneID][2].strip('"')
        ecod_x = ecodDB[geneID][3].strip('"')
        ecod_h = ecodDB[geneID][4].strip('"')
        ecod_t = ecodDB[geneID][5].strip('"')
        ecod_f = ecodDB[geneID][6].strip('"')
    else:
        ecod_id = "NaN"
        ecod_tm = "NaN"
        ecod_a = "NaN"
        ecod_x = "NaN"
        ecod_h = "NaN"
        ecod_t = "NaN"
        ecod_f = "NaN"
        
    if geneID in pdbDB:
        pdb_id = pdbDB[geneID][0].strip('"')
        pdb_tm = pdbDB[geneID][1].strip('"')
        pdb_descr = pdbDB[geneID][2].strip('"')
        pdb_classif = pdbDB[geneID][3].strip('"')
    else:
        pdb_id = "NaN"
        pdb_tm = "NaN"
        pdb_descr= "NaN"
        pdb_classif = "NaN"

    if geneID in af2_swissDB:
        af2_swiss_id = af2_swissDB[geneID][0].strip('"')
        af2_swiss_tm = af2_swissDB[geneID][1].strip('"')
        af2_swiss_name = af2_swissDB[geneID][2].strip('"')
        af2_swiss_function = af2_swissDB[geneID][3].strip('"')
    else:
        af2_swiss_id = "NaN"
        af2_swiss_tm = "NaN"
        af2_swiss_name = "NaN"
        af2_swiss_function = "NaN"

    if geneID in af2_proteomeDB:
        af2_prot_id = af2_proteomeDB[geneID][0].strip('"')
        af2_prot_tm = af2_proteomeDB[geneID][1].strip('"')
        af2_prot_name = af2_proteomeDB[geneID][2].strip('"')
        af2_prot_function = af2_proteomeDB[geneID][3].strip('"')
    else:
        af2_prot_id = "NaN"
        af2_prot_tm = "NaN"
        af2_prot_name = "NaN"
        af2_prot_function = "NaN"

    if geneID in dpi04_expressionDB:
        dpi04_l2fc = dpi04_expressionDB[geneID][0].strip('"')
        dpi04_padj = dpi04_expressionDB[geneID][1].strip('"')
    else:
        dpi04_l2fc = "NaN"
        dpi04_padj = "NaN"        

    if geneID in dpi07_expressionDB:
        dpi07_l2fc = dpi07_expressionDB[geneID][0].strip('"')
        dpi07_padj = dpi07_expressionDB[geneID][1].strip('"')
    else:
        dpi07_l2fc = "NaN"
        dpi07_padj = "NaN"

    if geneID in dpi10_expressionDB:
        dpi10_l2fc = dpi10_expressionDB[geneID][0].strip('"')
        dpi10_padj = dpi10_expressionDB[geneID][1].strip('"')
    else:
        dpi10_l2fc = "NaN"
        dpi10_padj = "NaN"


    # create a list with the values for this row
    row = [geneID, length, mature_length, secretome, disorder_perc, disorder_num, pLDDT_esmfold, prob_AMP_esmfold, annot_AMP_esmfold, pLDDT_af2, prob_AMP_af2, annot_AMP_af2, tmscore, rmsd, sp_start, sp_end, signalp, sp2, tm2, tmhmm, interpro, targetp, deeploc, effectorp, cazymes, proteases, lipases, bsc_id, bsc_function, bsc_gene_function, orthogroup, kog, go, pathway, desctrip, subgraph, mcode, cath_id, cath_tm, cath_class, cath_archit, cath_topo, cath_family, scope_id, scope_tm, scope_class, scope_fold, scope_superfam, scope_fam, ecod_id, ecod_tm, ecod_a, ecod_x, ecod_h, ecod_t, ecod_f,pdb_id, pdb_tm, pdb_descr, pdb_classif, af2_swiss_id, af2_swiss_tm, af2_swiss_name, af2_swiss_function, af2_prot_id, af2_prot_tm,  af2_prot_name, af2_prot_function, dpi04_l2fc, dpi04_padj, dpi07_l2fc, dpi07_padj, dpi10_l2fc, dpi10_padj, aa_seq]
   
    # add this row to the list of rows
    rows.append(row)

# write the table 
with out_anno as outfile:
    writer = csv.writer(outfile, delimiter='\t')
    writer.writerow(header)
    for row in rows:
        writer.writerow([str(i).replace("[","").replace("]","").replace("'","") if type(i)==list else i for i in row])





