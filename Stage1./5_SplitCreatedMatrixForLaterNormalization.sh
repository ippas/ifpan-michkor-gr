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

input_file="output_file4_withoutEID.csv"
output_directory="CHUNKdirectory"
chunk_size=500

# Create output directory if it doesn't exist
mkdir -p $output_directory

# Get the total number of columns in the input file
num_columns=$(awk -F, 'NR==1 {print NF}' $input_file)

# Loop through chunks
for ((start=1; start<=$num_columns; start+=chunk_size)); do
    end=$((start+chunk_size-1))
    
    # Create output file name based on chunk range
    output_file="${output_directory}/chunk_${start}_${end}.csv"

    # Extract columns using awk
    awk -v start=$start -v end=$end -F, '{
        for (i=start; i<=end; i++) {
            printf "%s%s", $i, (i<end ? "," : "");
        }
        printf "\n";
    }' $input_file > $output_file

    echo "Created chunk: $output_file"
done
exit;