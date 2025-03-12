#!/bin/bash

#BATCH --job-name=targetp2 #Give your job a name.

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


# Set the directory containing your protein files
protein_dir="/home/dalsasso/data/References/secretomes/"

# Set the output directory for Phobius results
targetp2_output_dir="/home/dalsasso/annotations/targetp2/"

# Iterate over protein files in the directory
for protein_file in "$protein_dir"*secretome.fa; do
   
    # Get the base name (before the first ".") of the file
    base_name=$(basename "$protein_file" | cut -d. -f1)

    # Run Targetp2
    /home/dalsasso/softwares/targetp-2.0/bin/targetp -gff3 -format short -org non-pl -fasta "$protein_file" -prefix "${targetp_output_dir}${base_name}.secretome.targetp2"

done
