import csv

# Read PhenotypesCodes.csv and create a dictionary mapping numbers to codes
number_to_code = {}
with open('combined_output.csv', 'r') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    for row in reader:
        for index, code in enumerate(row):
            number_to_code[index] = code.strip()

# Read FinalPhenoICD10_output.tsv and create a dictionary mapping codes to descriptions
code_to_description = {}
with open('FinalPhenoGWASandICD10_output.tsv', 'r') as tsvfile:
    reader = csv.reader(tsvfile, delimiter='\t')
    for row in reader:
        if len(row) >= 2:
            code_to_description[row[0]] = row[1]

# Read zlist.txt and extract corresponding codes and descriptions
########### zlist.txt prepared manually ################
output_lines = []
current_component = None
with open('zlistFreq50.txt', 'r') as input_file:
    lines = input_file.readlines()
    for line in lines:
        line = line.strip()
        if line.startswith('Component'):
            current_component = line
            output_lines.append(line)
        elif line.isdigit() and current_component is not None:
            number = int(line)
            if number in number_to_code:
                code = number_to_code[number]
                description = code_to_description.get(code, '')
                output_lines.append(f'{number}\t{code}\t{description}')
        else:
            output_lines.append(line)

# Write the output to zlist2_output.txt
with open('zlistFreq50_output.txt', 'w') as output_file:
    for line in output_lines:
        output_file.write(line + '\n')
