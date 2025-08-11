#!/bin/bash
#SBATCH --job-name=amapec
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=50G
#SBATCH --time=48:00:00
#SBATCH --output=job.%J.out
#SBATCH --error=job.%J.err

# Load modules:
module load gcc12-env/12.3.0
module load miniconda3/23.5.2
module load R/4.3.1 gcc/12.3.0


conda activate amapec_env


/gxfs_home/cau/sunbo511/amapec/amapec -i "/gxfs_work/cau/sunbo511/alphafold/Zpa796_AF2bestmodels/" -o "./" -d -t 20


conda deactivate
jobinfo
