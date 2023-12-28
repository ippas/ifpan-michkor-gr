#!/bin/bash

# Read the input file
input_file="filtered_statistics50.tsv"

# Extract disease codes and save them in a temporary file
awk '{print $1}' "$input_file" | tail -n +2 > temp_disease_codes.txt

# Convert disease codes to the desired format
awk '{printf "%s,", $1}' temp_disease_codes.txt | sed 's/,$/\n/' > output.txt

# Save the output to a CSV file
echo "$(cat output.txt)" > PhenotypesCodesFreq50_copy_no_eid.csv

# Clean up temporary files
rm temp_disease_codes.txt output.txt

echo "Output saved to PhenotypesCodesFreq50_copy_no_eid.csv"

head -1 outputNomissing.csv > header_output.txt

awk -F',' '{for(i=2;i<=33;i++) printf "%s%s", $i, (i<33?",":"\n")}' header_output.txt > columns_2_to_33.txt

mv columns_2_to_33.txt columns_2_to_33.csv

paste -d ',' columns_2_to_33.csv PhenotypesCodesFreq50_copy_no_eid.csv > combined_output.csv

rm columns_2_to_33.csv