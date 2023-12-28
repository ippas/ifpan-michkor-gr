import csv

# Define file paths
matches_file_path = "output_matches.txt"
input_csv_path = "FinalICD10Normalized.csv"
output_csv_path = "FINALofFINALsICD10Normalized77_2.csv"

# Read indices from the matches file
with open(matches_file_path, "r") as matches_file:
    match_indices = [int(line.strip()) for line in matches_file]

# Read the input CSV and extract matching rows
with open(input_csv_path, "r") as input_csv, open(output_csv_path, "w") as output_csv:
    for row_index, row in enumerate(input_csv, start=1):
        if row_index in match_indices:
            output_csv.write(row)

print(f"Matching rows extracted and saved to {output_csv_path}")
