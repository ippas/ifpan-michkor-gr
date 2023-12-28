import csv

# Define the line numbers to extract
# Row numbers to be +2 to match with correct phenotype 
line_numbers_to_extract = [
    6, 49, 55, 136, 152  # Updated line numbers
]

# Specify the input CSV file
csv_file_path = "Transposed_by_conda_transpose50.csv"
output_file_path = "extracted_lines_output.txt"

# Open the CSV file and write the extracted lines to a new file
with open(csv_file_path, 'r', newline='') as file:
    csv_reader = csv.reader(file)
    with open(output_file_path, 'w') as output_file:
        for i, row in enumerate(csv_reader, start=1):
            if i in line_numbers_to_extract:
                output_file.write(','.join(row) + '\n')

print(f"Extracted lines have been saved to {output_file_path}.")
