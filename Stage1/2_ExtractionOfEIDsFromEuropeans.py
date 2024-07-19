input_file = "NewData.txt"
output_file = "Europeans.txt"
target_values = {'1', '1001', '1002'}

# Read data from the input file
with open(input_file, 'r') as f:
    lines = f.readlines()

# Extract eids with target values along with their related data
matching_lines = []
for line in lines[1:]:  # Skip the header line
    columns = line.split()
    eid = columns[0]
    values = set(columns[1:])
    if any(value in target_values for value in values):
        matching_lines.append(line)

# Save the extracted lines to the output file
with open(output_file, 'w') as f:
    f.write("eid    3-0.0    21000-0.0    21000-1.0    21000-2.0    21000-3.0\n")
    for matching_line in matching_lines:
        f.write(matching_line)

print(f"Extracted data for matching eids saved to {output_file}")
