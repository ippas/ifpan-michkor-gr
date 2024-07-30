# import modules
import os
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from tqdm import tqdm
import multiprocessing as mp
import csv
from concurrent.futures import ProcessPoolExecutor


print(os.getcwd())

########################################################################################
# functions
########################################################################################
def standardize_dataframe(df):
    """
    Standardizes each column of the given DataFrame.
    
    Parameters:
    - df: A pandas DataFrame where each numerical column needs to be standardized.
    
    Returns:
    - A pandas DataFrame with each column standardized.
    """
    scaler = StandardScaler()
    for column in df.columns:
        df[column] = scaler.fit_transform(df[[column]])
    return df


def standardize_dataframe(df):
    """
    Standardizes each column of the given DataFrame, with a progress bar.
    
    Parameters:
    - df: A pandas DataFrame where each numerical column needs to be standardized.
    
    Returns:
    - A pandas DataFrame with each column standardized.
    """
    scaler = StandardScaler()
    for column in tqdm(df.columns, desc="Standardizing Columns"):
        df[column] = scaler.fit_transform(df[[column]])
    return df


def save_large_dataframe_to_tsv(df, output_filename, num_chunks=None):
    def save_chunk_to_tsv(chunk, chunk_index, base_filename):
        chunk_filename = f"{base_filename}_part_{chunk_index}.tsv"
        chunk.to_csv(chunk_filename, sep='\t', index=False)

    def combine_tsv_files(num_chunks, base_filename, final_filename):
        with open(final_filename, 'w') as outfile:
            for i in range(num_chunks):
                chunk_filename = f"{base_filename}_part_{i}.tsv"
                with open(chunk_filename, 'r') as infile:
                    if i != 0:
                        infile.readline()  # Skip header for all but the first file
                    outfile.write(infile.read())
                os.remove(chunk_filename)  # Remove chunk file after combining

    if num_chunks is None:
        num_chunks = mp.cpu_count()  # Use all available cores if num_chunks is not specified

    chunk_size = int(np.ceil(len(df) / num_chunks))
    chunks = [df.iloc[i*chunk_size:(i+1)*chunk_size] for i in range(num_chunks)]
    base_filename = "temp_chunk_output"

    with mp.Pool(processes=num_chunks) as pool:
        pool.starmap(save_chunk_to_tsv, [(chunk, i, base_filename) for i, chunk in enumerate(chunks)])

    combine_tsv_files(num_chunks, base_filename, output_filename)

########################################################################################
# Read the RealDiseaseCodesForPatients.txt file
########################################################################################

# read the RealDiseaseCodesForPatients.txt file
ukb_phenotypes_df = pd.read_csv('/net/pr2/projects/plgrid/plggneuromol/matzieb/projects/ifpan-michkor-gr/data/UKB-phenotypes/ukb-phenotypes-merged.tsv', sep='\t')

# from uk_phenotypes_df, select columns with 41270 prefix in column names
columns_41270 = [col for col in ukb_phenotypes_df.columns if '41270' in col]
print(columns_41270)

# select columns with 41270 prefix and eid column
columns_41270_eid = ['eid'] + columns_41270
ukb_phenotypes_41270_df = ukb_phenotypes_df[columns_41270_eid]

# save the selected columns to a file
ukb_phenotypes_41270_df.to_csv('/net/pr2/projects/plgrid/plggneuromol/matzieb/projects/ifpan-michkor-gr/data/UKB-phenotypes/ukb-phenotypes-41270-icd10.tsv', sep='\t', index=False)

########################################################################################
# Read the ICD10 codes
########################################################################################

# Load the data specifying 'eid' as the index column
file_path = '/net/pr2/projects/plgrid/plggneuromol/matzieb/projects/ifpan-michkor-gr/data/UKB-phenotypes/ukb-phenotypes-41270-icd10.tsv'
ukb_icd10_41270_df = pd.read_csv(file_path, sep='\t', index_col='eid')

# Rename the columns for easier handling
ukb_icd10_41270_df.columns = ukb_icd10_41270_df.columns.str.replace('41270_', '')

# Flatten the DataFrame to get all unique ICD10 codes, excluding NaN values
unique_icd10_codes = pd.Series(ukb_icd10_41270_df.values.ravel()).dropna().unique()

# Create an index mapping for each unique code
icd10_index = {code: idx for idx, code in enumerate(unique_icd10_codes)}

# Prepare an empty matrix for one-hot encoding
icd10_one_hot_matrix = np.zeros((len(ukb_icd10_41270_df), len(unique_icd10_codes)), dtype=int)

# Populate the matrix
for idx, (eid, codes) in enumerate(ukb_icd10_41270_df.iterrows()):
    valid_codes = codes.dropna()
    indices = [icd10_index[code] for code in valid_codes if code in icd10_index]
    icd10_one_hot_matrix[idx, indices] = 1

# Convert the matrix to a DataFrame for better manipulation and visualization
icd10_one_hot_df = pd.DataFrame(icd10_one_hot_matrix, index=ukb_icd10_41270_df.index, columns=unique_icd10_codes)

icd10_one_hot_df.iloc[1:10, 1:10]

icd10_one_hot_df["E780"].tolist()

# Save the one-hot encoded DataFrame to a file
# icd10_one_hot_df.to_csv('/net/pr2/projects/plgrid/plggneuromol/matzieb/projects/ifpan-michkor-gr/data/UKB-phenotypes/ukb-phenotypes-41270-icd10-one-hot.tsv', sep='\t')




########################################################################################
# standrdize the ICD10 codes
########################################################################################

icd10_one_hot_standardize_df = standardize_dataframe(icd10_one_hot_df)

icd10_one_hot_standardize_df["E785"].shape

icd10_one_hot_standardize_df.tail()

# move index to a column 
icd10_one_hot_standardize_df.reset_index(inplace=True)

# eid as character
icd10_one_hot_standardize_df['eid'] = icd10_one_hot_standardize_df['eid'].astype(str)

# detect available cores
mp.cpu_count()


########################################################################################
# save the standardized ICD10 codes. Test 1.
def save_chunk_to_tsv(chunk, filename, mode):
    """
    Save a chunk of a DataFrame to a TSV file using the csv module.
    
    Parameters:
    chunk (pd.DataFrame): The chunk of the DataFrame to save.
    filename (str): The name of the TSV file to save the chunk to.
    mode (str): The file mode ('w' for write, 'a' for append).
    """
    with open(filename, mode=mode, newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        if mode == 'w':
            writer.writerow(chunk.columns)  # Write header
        writer.writerows(chunk.values)

def save_df_to_tsv_custom(df, filename, chunksize=50000, n_cores=4):
    """
    Save a DataFrame to a TSV file using the csv module for better performance, with a progress bar and multi-core processing.
    
    Parameters:
    df (pandas.DataFrame): The DataFrame to save.
    filename (str): The name of the TSV file to save the DataFrame to.
    chunksize (int): Number of rows per chunk to write at a time.
    n_cores (int): Number of cores to use for parallel processing.
    """
    total = len(df)
    chunks = [df.iloc[i:i + chunksize] for i in range(0, total, chunksize)]
    
    with tqdm(total=total, desc="Saving DataFrame to TSV", unit='rows') as pbar:
        with ProcessPoolExecutor(max_workers=n_cores) as executor:
            futures = []
            for i, chunk in enumerate(chunks):
                mode = 'w' if i == 0 else 'a'
                futures.append(executor.submit(save_chunk_to_tsv, chunk, filename, mode))
                pbar.update(len(chunk))
            
            for future in futures:
                future.result()  # Ensure all futures have completed

# save_df_to_tsv_custom(icd10_one_hot_stnadardize_df, 
#                       "data/UKB-phenotypes/icd10-41270-one-hot-standardized.tsv", 
#                       chunksize=5000, n_cores=24)

########################################################################################
# save the standardized ICD10 codes. Test 2.

import pandas as pd
import csv
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor

def save_chunk_to_tsv(chunk, filename, mode):
    """
    Save a chunk of a DataFrame to a TSV file using the csv module.
    
    Parameters:
    chunk (pd.DataFrame): The chunk of the DataFrame to save.
    filename (str): The name of the TSV file to save the chunk to.
    mode (str): The file mode ('w' for write, 'a' for append).
    """
    with open(filename, mode=mode, newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        if mode == 'w':
            writer.writerow(chunk.columns)  # Write header
        writer.writerows(chunk.values)

def save_df_to_tsv_custom(df, filename, chunksize=50000, n_cores=4):
    """
    Save a DataFrame to a TSV file using the csv module for better performance, with a progress bar and multi-core processing.
    
    Parameters:
    df (pandas.DataFrame): The DataFrame to save.
    filename (str): The name of the TSV file to save the DataFrame to.
    chunksize (int): Number of rows per chunk to write at a time.
    n_cores (int): Number of cores to use for parallel processing.
    """
    # Add an index column to preserve the original order
    df['index'] = df.index
    
    total = len(df)
    chunks = [df.iloc[i:i + chunksize] for i in range(0, total, chunksize)]
    
    with tqdm(total=total, desc="Saving DataFrame to TSV", unit='rows') as pbar:
        with ProcessPoolExecutor(max_workers=n_cores) as executor:
            futures = []
            for i, chunk in enumerate(chunks):
                mode = 'w' if i == 0 else 'a'
                futures.append(executor.submit(save_chunk_to_tsv, chunk.sort_values(by='index'), filename, mode))
                pbar.update(len(chunk))
            
            for future in futures:
                future.result()  # Ensure all futures have completed

    # Optionally, remove the index column from the final file
    df.drop(columns=['index'], inplace=True)

# Example usage
# save_df_to_tsv_custom(icd10_one_hot_stnadardize_df, 
#                       "data/UKB-phenotypes/icd10-41270-one-hot-standardized2.tsv", 
#                       chunksize=25000, n_cores=10)

########################################################################################
# save the standardized ICD10 codes. Test 3.
import pandas as pd
import csv
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor

def save_chunk_to_tsv(chunk, filename, mode):
    """
    Save a chunk of a DataFrame to a TSV file using the csv module.
    
    Parameters:
    chunk (pd.DataFrame): The chunk of the DataFrame to save.
    filename (str): The name of the TSV file to save the chunk to.
    mode (str): The file mode ('w' for write, 'a' for append).
    """
    with open(filename, mode=mode, newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        if mode == 'w':
            writer.writerow(chunk.columns)  # Write header
        writer.writerows(chunk.values)

def save_df_to_tsv_custom(df, filename, chunksize=50000, n_cores=4):
    """
    Save a DataFrame to a TSV file using the csv module for better performance, with a progress bar and multi-core processing.
    
    Parameters:
    df (pandas.DataFrame): The DataFrame to save.
    filename (str): The name of the TSV file to save the DataFrame to.
    chunksize (int): Number of rows per chunk to write at a time.
    n_cores (int): Number of cores to use for parallel processing.
    """
    # Add an index column to preserve the original order
    df.reset_index(drop=False, inplace=True)
    
    total = len(df)
    mode = 'w'
    
    with tqdm(total=total, desc="Saving DataFrame to TSV", unit='rows') as pbar:
        for i in range(0, total, chunksize):
            chunk = df.iloc[i:i + chunksize]
            save_chunk_to_tsv(chunk, filename, mode)
            mode = 'a'  # Change mode to append after the first chunk
            pbar.update(len(chunk))

    # Optionally, remove the index column from the DataFrame
    df.drop(columns=['index'], inplace=True)
    
save_df_to_tsv_custom(icd10_one_hot_stnadardize_df, 
                      "data/UKB-phenotypes/icd10-41270-one-hot-standardized.tsv", 
                      chunksize=5000, n_cores=24)

# ########################################################################################
# # save the standardized ICD10 codes. Test 4.

# import dask.dataframe as dd
# from dask.diagnostics import ProgressBar
# import pandas as pd


# def save_and_merge_to_tsv(df, output_dir='./', npartitions=10, combined_filename='output_combined.tsv'):
#     """
#     Save a pandas DataFrame to multiple TSV files using Dask with progress bar and merge them into one file.

#     Parameters:
#     - df: pandas DataFrame to be saved.
#     - output_dir: Directory to save the TSV files.
#     - npartitions: Number of partitions to split the DataFrame for parallel processing.
#     - combined_filename: Name of the final combined TSV file.
#     """
#     # Convert pandas DataFrame to Dask DataFrame
#     ddf = dd.from_pandas(df, npartitions=npartitions)

#     # Save the Dask DataFrame to TSV files with a progress bar
#     with ProgressBar():
#         ddf.to_csv(os.path.join(output_dir, 'output-*.tsv'), sep='\t', index=False)

#     # Get list of all TSV files in the directory
#     tsv_files = sorted([f for f in os.listdir(output_dir) if f.startswith('output-') and f.endswith('.tsv')])

#     # Combine all TSV files into one
#     combined_file_path = os.path.join(output_dir, combined_filename)
#     with open(combined_file_path, 'w') as outfile:
#         for i, fname in enumerate(tsv_files):
#             with open(os.path.join(output_dir, fname)) as infile:
#                 if i != 0:
#                     # Skip the header in subsequent files
#                     next(infile)
#                 outfile.write(infile.read())

#     print(f"All TSV files have been merged into {combined_file_path}.")
    
    

save_and_merge_to_tsv(icd10_one_hot_stnadardize_df.head(50000),
                      output_dir='data/',
                      npartitions=20,
                      combined_filename='icd10-41270-one-hot-standardized3.tsv')
