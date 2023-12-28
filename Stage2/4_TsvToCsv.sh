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
# Extract lines
# sed -n '150001,300000p' New_Data.tsv > 150001_300000.tsv

# Convert tsv to csv
## sed 's/\t/,/g' NonMissingData2.tsv > NonMissingData2.csv

# Fill empty cells in original data with 0
# awk 'BEGIN {FS=OFS=","} {for(i=1; i<=NF; i++) if($i=="") $i=0; print}' NonMissingData2.csv > outputNomissing.csv
## awk 'BEGIN {FS=OFS=","} {for(i=1; i<=NF; i++) if($i=="") $i=0; print}' Test1.csv > outputNomissing2.csv
# Replace NA with 0
# NA may raised during data standardization
## sed -i 's/NA/0/g' outputNomissing2.csv

# Removing first row and first column
tail -n +2 Transposed_by_conda_transpose500.csv > New_Normalized_Data.csv
cut -d',' -f2- New_Normalized_Data.csv > New_DataFreq500.csv
rm New_Normalized_Data.csv

#### Making for loop for all Test*.csv files ####
# Loop through all CSV files starting with "Test"
# for file in Test*.csv; do
#     # Extract file name without extension
#     filename=$(basename -- "$file")
#     filename_noext="${filename%.*}"
#     
#     # Remove header line and create new CSV file
#     tail -n +2 "$file" > "New_$filename_noext.csv"
#     
#     # Optional: Print a message for each processed file
#     echo "Processed $file"
# done
# 
# echo "All CSV files processed successfully."

exit;