#!/bin/bash

#BATCH --job-name=interproscan #Give your job a name.

#SBATCH --nodes=1 #Only increase for openmpi jobs.
#SBATCH --ntasks=1 #e.g if set to 2, could run two softwares in the script at the same time.
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6 #Multithreading.
#SBATCH --time=72:00:00 #Time for a job to run given as hh:mm:ss.
#SBATCH --mem=10G #Total Memory per node to use for the job
#SBATCH --error=job.%J.err #Std Error write standard error to this file
#SBATCH --output=job.%J.out #Std Out write standard output to this file
#SBATCH --mail-type=ALL #Notify user by email when certain event types occur (BEGIN, END, FAIL, REQUEUE)
#SBATCH --mail-user=dalsasso@evolbio.mpg.de #Email for notifications from previous line
#SBATCH --partition=standard #Request a specific partition for the resource allocation.


##########################
# Module load:
module load java/x64/17u3

#########################
# Run commands from here:

# Version v5.54-87.0

protein_dir="/home/dalsasso/data/References/proteomes/"
interpro_output_dir="/home/dalsasso/annotations/interpro/"


for protein_file in "$protein_dir"*no_stop_codon.fa; do
   
    base_name=$(basename "$protein_file" | cut -d_ -f1)

    interproscan -cpu 6 -appl SMART,SUPERFAMILY,CDD,TIGRFAM,Pfam,Coils,Gene3D \
         -i "$protein_file" -b "${interpro_output_dir}${base_name}.interpro" -f tsv
done





