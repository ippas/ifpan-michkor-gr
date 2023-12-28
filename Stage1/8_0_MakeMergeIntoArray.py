import csv

# Initialize an empty list to store the merged data
merged_data = []

# Columns to merge into arrays
columns_to_merge = ['2453-0.0', '2453-1.0', '2453-2.0', '2453-3.0']

# Path to TSV file
tsv_file_path = '2453.tsv'

# Open the TSV file and read its contents
with open(tsv_file_path, newline='') as tsvfile:
    reader = csv.DictReader(tsvfile, delimiter='\t')
    
    # Iterate through each row in the TSV file
    for row in reader:
        # Initialize a dictionary for the merged row
        merged_row = {}
        
        # Iterate through columns in the row
        for column, value in row.items():
            # If the column is in the columns_to_merge list, split and store as an array
            if column in columns_to_merge:
                merged_row[column] = [cell.strip() for cell in value.split(',')]
            else:
                merged_row[column] = value.strip()
        
        # Append the merged row to the merged_data list
        merged_data.append(merged_row)

# Print the merged data as dictionaries
for row in merged_data:
    print(row)

# Define a new file path for the output TSV file
output_file_path = '2453_merged_data.tsv'

# Extract the column names from the first row of the input file
column_names = merged_data[0].keys()

# Open the output TSV file for writing
with open(output_file_path, 'w', newline='') as outfile:
    writer = csv.DictWriter(outfile, delimiter='\t', fieldnames=column_names)
    
    # Write the merged data to the output TSV file
    writer.writeheader()
    writer.writerows(merged_data)

# Print a message to confirm that the data has been written to the output file
print(f"Merged data has been written to {output_file_path}")
