#!/bin/bash
###############################################################
### SLURM preamble:  ###
###############################################################
#SBATCH --account plgwgsdepresja3-cpu
#SBATCH --partition plgrid
#SBATCH --nodes=1
#SBATCH --mem=180G
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time 72:00:00
#SBATCH -C localfs
#############################################################

# head -300000 outputNomissingFreq500.csv > outputNomissingFreq500_For300000lines.csv

# Convert csv to tsv
sed 's/,/\t/g' ReadyForPCAFreq50.csv > ReadyForPCAFreq50.tsv

# Transpose the tsv file
source activate Transpose
transpose -s <ReadyForPCAFreq50.tsv > Transposed_by_conda_transpose50.tsv
conda deactivate

# Convert tsv to csv
sed 's/\t/,/g' Transposed_by_conda_transpose50.tsv > Transposed_by_conda_transpose50.csv


exit;