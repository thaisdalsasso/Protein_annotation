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


#### Set the directory containing your protein files
###protein_dir="/Users/dalsasso/Desktop/Posdoc/CAU/data/References/secretomes/"
###
#### Set the output directory for Phobius results
###deeploc2_output_dir="/Users/dalsasso/Desktop/Posdoc/CAU/annotations/deeploc2"
###
#### Iterate over protein files in the directory
###for protein_file in "$protein_dir"*secretome.fa; do
###   
###    # Get the base name (before the first ".") of the file
###    base_name=$(basename "$protein_file" | cut -d. -f1)
###
###    # Run DeepLoc2
###    /Users/dalsasso/Library/Python/3.9/bin/deeploc2 -f "$protein_file" -m Fast -o "${deeploc2_output_dir}${base_name}.secretome.deeploc2"
###
###done


/Users/dalsasso/Library/Python/3.9/bin/deeploc2 -f ../../data/References/secretomes/CbeticolaCb09-40.secretome.fa -m Fast -o CbeticolaCb09-40
/Users/dalsasso/Library/Python/3.9/bin/deeploc2 -f ../../data/References/secretomes/Ndiscreta.secretome.fa -m Fast -o Ndiscreta 
/Users/dalsasso/Library/Python/3.9/bin/deeploc2 -f ../../data/References/secretomes/Pteresteres0-1.secretome.fa -m Fast -o Pteresteres0-1
/Users/dalsasso/Library/Python/3.9/bin/deeploc2 -f ../../data/References/secretomes/Zpa796.secretome.fa -m Fast -o Zpa796
/Users/dalsasso/Library/Python/3.9/bin/deeploc2 -f ../../data/References/secretomes/Zt469.secretome.fa -m Fast -o Zt469  