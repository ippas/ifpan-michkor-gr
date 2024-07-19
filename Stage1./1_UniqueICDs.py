import pandas as pd

# Read the original TSV file into a DataFrame
df = pd.read_csv('ICD_10s.tsv', sep='\t')

# Remove duplicates from the 'DiseaseCode' column and create a new DataFrame
unique_df = df.drop_duplicates(subset='DiseaseCode')

# Save the unique DataFrame to a new TSV file
unique_df.to_csv('unique_disease_codes.tsv', sep='\t', index=False)
