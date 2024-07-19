import pandas as pd

# Load the TSV file into a DataFrame
df = pd.read_csv('2492_merged_data.tsv', sep='\t')

# Iterate through each column that I want to one-hot encode
columns_to_encode = ['2492-0.0', '2492-1.0', '2492-2.0', '2492-3.0']

for column in columns_to_encode:
    # Extract the unique values in the column
    unique_values = df[column].unique()
    
    # Iterate through unique values and create one-hot encoded columns
    for value in unique_values:
        encoded_column = f"{column}_{value}"
        df[encoded_column] = (df[column] == value).astype(int)
    
# Drop the original columns that were one-hot encoded
df.drop(columns=columns_to_encode, inplace=True)

# Save the one-hot encoded DataFrame to a new TSV file
df.to_csv('2492_one_hot_encoded_data.tsv', sep='\t', index=False)

# Print the one-hot encoded data
print(df)
