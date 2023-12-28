import numpy as np

# Read the disease codes for patients file
patient_data = np.genfromtxt("RealDiseaseCodesForPatients2.txt", delimiter='\t', dtype=str)
eids = patient_data[:, 0].astype(int)
disease_codes = patient_data[:, 1:]

# Read the unique disease codes file
unique_codes = np.genfromtxt("unique_disease_codes.tsv", delimiter='\t', dtype=str, skip_header=1)

# Create a mapping between unique disease codes and column indices in the output matrix
code_to_index = {code: i for i, code in enumerate(unique_codes)}

# Create a matrix with sample IDs as rows and unique disease codes as columns
output_matrix = np.column_stack((eids, np.zeros((len(eids), len(unique_codes)), dtype=int)))

# Iterate through samples and update the output matrix
for i, row_codes in enumerate(disease_codes):
    for code in row_codes:
        if code in code_to_index:
            output_matrix[i, code_to_index[code] + 1] = 1

# Save the final output to a new file
header = 'eid,' + ','.join(unique_codes)
np.savetxt("output_file4.csv", output_matrix, fmt='%d', delimiter=',', header=header, comments='')
