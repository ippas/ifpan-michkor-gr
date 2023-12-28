from collections import Counter

# Read the RealDiseaseCodesForPatients2.txt file
with open("RealDiseaseCodesForPatients2.txt", "r") as file:
    lines = file.readlines()

# Extract disease codes from each line and count their frequencies
all_disease_codes = [code.strip() for line in lines for code in line.split('\t')[1:] if code.strip()]
disease_code_counts = Counter(all_disease_codes)

# Read the unique_disease_codes.tsv file
with open("unique_disease_codes.tsv", "r") as file:
    unique_disease_codes = [line.strip() for line in file.readlines()]

# Filter and create a new file with unique disease codes and their frequencies
output_lines = []
for code in unique_disease_codes:
    frequency = disease_code_counts.get(code, 0)
    output_lines.append(f"{code}\t{frequency}")

# Write the output to a new file
with open("unique_disease_code_frequencies.txt", "w") as output_file:
    output_file.write("\n".join(output_lines))
