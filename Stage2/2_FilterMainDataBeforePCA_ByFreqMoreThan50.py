import pandas as pd

# Define the file paths
input_file_path = 'ReadyForPCA.csv'
filtered_stats_file_path = 'filtered_statistics500.tsv'

# Read the relevant columns from the large CSV file
columns_to_extract = ['eid', '31-0.0', '34-0.0', '48-', '49-', '137-', '4079-', '4080-', '12144-', '21001-', '21002-', '30000-', '30010-', '30040-', '30050-', '30080-', '30100-', '30130-', '40008-', '1259_', '20551_', '2492_', '6153_', '6156_', '1558_', '2443_', '4041_', '6154_', '6177_', '1618_', '2453_', '54_', '6155_']

# Additional columns from the filtered statistics file
filtered_stats_df = pd.read_csv(filtered_stats_file_path, sep='\t')
additional_columns = filtered_stats_df['DiseaseCode'].tolist()

# Combine all columns to extract
columns_to_extract += additional_columns

# Read the large CSV file with only the selected columns
df = pd.read_csv(input_file_path, usecols=columns_to_extract)

# Save the result to a new CSV file
output_file_path = 'ReadyForPCAFreq500.csv'
df.to_csv(output_file_path, index=False)

print(f"Extracted data saved to {output_file_path}")
