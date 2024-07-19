import pandas as pd
from sklearn.preprocessing import StandardScaler

# Specify all the file names
file_names = [
    'chunk_1_500.csv',
    'chunk_501_1000.csv',
    'chunk_1001_1500.csv',
    'chunk_1501_2000.csv',
    'chunk_2001_2500.csv',
    'chunk_2501_3000.csv',
    'chunk_3001_3500.csv',
    'chunk_3501_4000.csv',
    'chunk_4001_4500.csv',
    'chunk_4501_5000.csv',
    'chunk_5001_5500.csv',
    'chunk_5501_6000.csv',
    'chunk_6001_6500.csv',
    'chunk_6501_7000.csv',
    'chunk_7001_7500.csv',
    'chunk_7501_8000.csv',
    'chunk_8001_8500.csv',
    'chunk_8501_9000.csv',
    'chunk_9001_9500.csv',
    'chunk_9501_10000.csv',
    'chunk_10001_10500.csv',
    'chunk_10501_11000.csv',
    'chunk_11001_11500.csv',
    'chunk_11501_12000.csv',
    'chunk_12001_12500.csv',
    'chunk_12501_13000.csv',
    'chunk_13001_13500.csv',
    'chunk_13501_14000.csv',
]

# Iterate through each file
for file_name in file_names:
    # Read the CSV file into a pandas DataFrame
    df = pd.read_csv(file_name)

    # Extract and save the header row
    header = df.columns.tolist()

    # Standardize (normalize) the data, excluding the header row
    scaler = StandardScaler()
    df.iloc[:, :] = scaler.fit_transform(df.iloc[:, :])

    # Write the standardized data back to a new CSV file
    output_file = file_name.replace('.csv', '_standardized.csv')
    df.to_csv(output_file, index=False, header=header)

print("Data standardization complete.")
