#!/bin/bash
# Job name:
#SBATCH --job-name=dada
# Project:
#SBATCH --account=nn9320k

# Wall clock limit:
#SBATCH --time=04:55:00
# Number of tasks (MPI ranks):
#SBATCH --ntasks=16
# Max memory usage per task:
#SBATCH --mem-per-cpu=60G

## Set up job environment:
source /cluster/bin/jobsetup
module purge # clear any inherited modules
set -o errexit # exit on errors

module load R/3.4.4
Rscript dada2algor.R

