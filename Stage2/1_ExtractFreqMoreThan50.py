import pandas as pd

# Read the TSV file into a DataFrame
file_path = 'unique_disease_code_frequencies.tsv'
df = pd.read_csv(file_path, sep='\t')

# Filter rows where Freq is more than 50
filtered_df = df[df['Freq'] > 500]

# Save the filtered DataFrame to a new TSV file
output_file_path = 'filtered_statistics500.tsv'
filtered_df.to_csv(output_file_path, sep='\t', index=False)

print(f"Filtered data saved to {output_file_path}")
