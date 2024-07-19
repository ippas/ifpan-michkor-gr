import csv

# Read the first column of both CSV files
with open("AllGWASphenotypesNormalized.csv", "r") as file:
    all_gwas_data = set(row[0] for row in csv.reader(file))

with open("output_file4.csv", "r") as file:
    output_file4_data = list(csv.reader(file))

# Find exact matches and save the row numbers to an output text file
matches_index = []
for i, row in enumerate(output_file4_data):
    if row[0] in all_gwas_data:
        matches_index.append(i + 1)

with open("output_matches.txt", "w") as file:
    for index in matches_index:
        file.write(f"{index}\n")

print("Matching rows saved to output_matches.txt")
