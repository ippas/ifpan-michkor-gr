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

source activate ZSTD
zstd -k -19 New_DataFreq50.csv -o New_DataFreq50.zst
zstd -k -19 New_DataFreq500.csv -o New_DataFreq500.zst
conda deactivate

exit;