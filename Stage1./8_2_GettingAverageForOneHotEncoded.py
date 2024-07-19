# Initialize an empty list to store the averages
averages = []

# Read the TSV file
with open('20551_filtered_data.tsv', 'r') as file:
    lines = file.readlines()

# Iterate through each line starting from the second line
for line in lines[1:]:
    # Split the line by tab to get individual values
    values = line.strip().split('\t')
    
    # Convert the values to integers
    values = [int(value) for value in values]
    
    # Calculate the average of the values for this line
    average = sum(values) / len(values)
    
    # Append the average to the list
    averages.append(average)

# Write the averages to a new TSV file
with open('20551_result_averaged.tsv', 'w') as output_file:
    # Write the header line
    output_file.write('\t'.join(['20551_']) + '\n')
    
    # Write the exact averages for each line
    for average in averages:
        output_file.write(f'{average}\n')

print("Averages written to '20551_result_averaged.tsv'")
