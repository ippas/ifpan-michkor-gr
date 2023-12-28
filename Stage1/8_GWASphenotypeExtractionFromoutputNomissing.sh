#!/bin/bash
###############################################################
### SLURM preamble:  ###
###############################################################
#SBATCH --account plgwgsdepresja3-cpu
#SBATCH --partition plgrid
#SBATCH --nodes=1
#SBATCH --mem=150G
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time 72:00:00
#SBATCH -C localfs
#############################################################

cut -d ',' -f 1-33 outputNomissing.csv > AllGWASphenotypesNormalized.csv

exit;