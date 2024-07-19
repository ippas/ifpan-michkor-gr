#!/bin/bash
###############################################################
### SLURM preamble:  ###
###############################################################
#SBATCH --account plgwgsdepresja3-cpu
#SBATCH --partition plgrid
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time 24:00:00
#SBATCH -C localfs
#############################################################

cat output_ukb.txt | cut -f 2,15811-16069 > RealDiseaseCodesForPatients.txt

exit;