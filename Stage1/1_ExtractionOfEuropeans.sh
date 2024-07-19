#!/bin/bash
###############################################################
### SLURM preamble:  ###
###############################################################
#SBATCH --account plgwgsdepresja3-cpu
#SBATCH --partition plgrid-now
#SBATCH --nodes=1
#SBATCH --mem=50G
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time 12:00:00
#SBATCH -C localfs
#############################################################

cut -f 1,2,9744-9747 output_ukb.txt > NewData.txt

exit;