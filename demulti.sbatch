#!/bin/bash
# Job name:
#SBATCH --job-name=dada
# Project:
#SBATCH --account=nn9320k

# Wall clock limit:
#SBATCH --time=02:55:00
# Number of tasks (MPI ranks):
#SBATCH --ntasks=1
# Max memory usage per task:
#SBATCH --mem-per-cpu=5G

## Set up job environment:
source /cluster/bin/jobsetup
module purge # clear any inherited modules
set -o errexit # exit on errors

module load python3 #load python

python3 Demultiplexer.py

