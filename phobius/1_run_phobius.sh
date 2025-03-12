#!/bin/bash

#BATCH --job-name=phobius #Give your job a name.

#SBATCH --nodes=1 #Only increase for openmpi jobs.
#SBATCH --ntasks=1 #e.g if set to 2, could run two softwares in the script at the same time.
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10 #Multithreading.
#SBATCH --time=72:00:00 #Time for a job to run given as hh:mm:ss.
#SBATCH --mem=10G #Total Memory per node to use for the job
#SBATCH --error=job.%J.err #Std Error write standard error to this file
#SBATCH --output=job.%J.out #Std Out write standard output to this file
#SBATCH --mail-type=FAIL #Notify user by email when certain event types occur (BEGIN, END, FAIL, REQUEUE)
#SBATCH --mail-user=dalsasso@evolbio.mpg.de #Email for notifications from previous line
#SBATCH --partition=standard #Request a specific partition for the resource allocation.


##########################
# Module load:
#module load $MODULE_NAME 

#########################
# Run commands from here:

protein_dir="/home/dalsasso/data/References/proteomes/"
phobius_output_dir="/home/dalsasso/annotations/phobius/"


for protein_file in "$protein_dir"*no_stop_codon.fa; do
   
    base_name=$(basename "$protein_file" | cut -d_ -f1)

    perl /home/dalsasso/softwares/phobius/phobius.pl -short "$protein_file" > "${phobius_output_dir}${base_name}.phobius.short.out"

done


