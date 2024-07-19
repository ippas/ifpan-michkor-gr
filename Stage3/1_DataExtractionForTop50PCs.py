import pandas as pd

# Read the phenotype codes from top20PCs.txt
with open("top50PCs.txt", "r") as file:
    phenotype_codes = [line.strip() for line in file]

# Read the adjusted genetic data
genetic_data_df = pd.read_csv("150000adjusted_genetic_data.tsv", sep="\t")

# Filter the genetic data for the phenotype codes that are present in the column headers
valid_phenotype_codes = [code for code in phenotype_codes if code in genetic_data_df.columns]
filtered_data_df = genetic_data_df[valid_phenotype_codes]
filtered_data_df = filtered_data_df.loc[:,~filtered_data_df.columns.duplicated()]

# Save the filtered data to a new TSV file
filtered_data_df.to_csv("input_data50.tsv", sep="\t", index=False)

