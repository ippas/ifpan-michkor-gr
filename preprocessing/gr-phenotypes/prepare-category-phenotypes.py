# import modules
import os
import pandas as pd
import numpy as np
import requests
from bs4 import BeautifulSoup
import time
import multiprocessing as mp
import polars as pl
import tqdm

################################################################################
# fucntions
################################################################################
def read_file_in_chunks(file_path, sep='\t', chunk_size=10000, max_rows=None):
    """
    Read a large file in chunks and combine them into a single DataFrame.
    Optionally, limit the maximum number of rows to be read.
    Measure and print the execution time.

    Parameters:
    file_path (str): The path to the file to be read.
    sep (str): The delimiter for the file (default is '\t').
    chunk_size (int): The size of each chunk (default is 10000).
    max_rows (int or None): The maximum number of rows to read (default is None, which reads the entire file).

    Returns:
    DataFrame: The combined DataFrame up to the maximum number of rows.
    """
    start_time = time.time()  # Start the timer

    # Initialize an empty list to hold the chunks
    chunks = []
    rows_read = 0
    chunk_number = 0

    # Read the file in chunks
    for chunk in pd.read_csv(file_path, sep=sep, chunksize=chunk_size):
        chunks.append(chunk)
        rows_read += len(chunk)
        chunk_number += 1
        chunk_time = time.time() - start_time
        print(f"Reading chunk {chunk_number}: rows {rows_read - len(chunk)} to {rows_read}, time elapsed: {chunk_time:.2f} seconds")

        if max_rows is not None and rows_read >= max_rows:
            break

    # Combine the chunks into a single DataFrame
    combined_df = pd.concat(chunks, axis=0)

    # Trim the combined DataFrame to the maximum number of rows if necessary
    if max_rows is not None and len(combined_df) > max_rows:
        combined_df = combined_df.iloc[:max_rows]

    end_time = time.time()  # End the timer
    execution_time = end_time - start_time  # Calculate the total execution time

    print(f"Total execution time: {execution_time:.2f} seconds")

    return combined_df

    """
    Read a large file in chunks using parallel processing and combine them into a single DataFrame.
    Optionally, limit the maximum number of rows to be read.
    Measure and print the execution time.

    Parameters:
    file_path (str): The path to the file to be read.
    sep (str): The delimiter for the file (default is '\t').
    chunk_size (int): The size of each chunk (default is 10000).
    max_rows (int or None): The maximum number of rows to read (default is None, which reads the entire file).
    num_cores (int): The number of CPU cores to use for parallel processing (default is 1).

    Returns:
    DataFrame: The combined DataFrame up to the maximum number of rows.
    """
    start_time = time.time()  # Start the timer

    # Determine the number of rows in the file
    total_rows = sum(1 for _ in open(file_path)) - 1  # Subtract 1 for the header row
    if max_rows is not None:
        total_rows = min(total_rows, max_rows)

    # Create a list of starting rows for each chunk
    start_rows = list(range(1, total_rows, chunk_size))  # Start from 1 to skip header

    # Create a pool of workers with the specified number of cores
    pool = mp.Pool(num_cores)

    # Read the chunks in parallel
    chunks = pool.starmap(read_chunk, [(file_path, sep, chunk_size, start_row) for start_row in start_rows])

    # Close the pool and wait for the work to finish
    pool.close()
    pool.join()

    # Sort the chunks by their starting row to ensure correct order
    chunks.sort(key=lambda x: x[0])

    # Extract the DataFrames from the sorted list
    ordered_chunks = [chunk for _, chunk in chunks]

    # Combine the chunks into a single DataFrame
    combined_df = pd.concat(ordered_chunks, axis=0)

    # Trim the combined DataFrame to the maximum number of rows if necessary
    if max_rows is not None and len(combined_df) > max_rows:
        combined_df = combined_df.iloc[:max_rows]

    end_time = time.time()  # End the timer
    execution_time = end_time - start_time  # Calculate the total execution time

    print(f"Total execution time: {execution_time:.2f} seconds")

    return combined_df

def get_biobank_metadata(field_id, number_attempts=3, time_duration=5):
    url = f"https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id={field_id}"
    
    for attempt in range(number_attempts):
        try:
            response = requests.get(url)
            response.raise_for_status()  # Raise an error for bad status codes
            soup = BeautifulSoup(response.content, 'html.parser')

            # Extract description
            description_row = soup.find('td', text="Description:")
            description = description_row.find_next_sibling('td').text.strip() if description_row else None

            # Extract category
            category_row = soup.find('td', text="Category:")
            category = []
            if category_row:
                category = [a.text for a in category_row.find_next_sibling('td').find_all('a')]

            metadata = {
                "Description": description,
                "Category": category,
                "Url": url,
            }

            return metadata
        
        except requests.RequestException as e:
            print(f"Attempt {attempt + 1} failed: {e}")
            if attempt < number_attempts - 1:
                time.sleep(time_duration)
            else:
                raise
            
def get_biobank_metadata(field_id, number_attempts=3, time_duration=5):
    url = f"https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id={field_id}"
    
    for attempt in range(number_attempts):
        try:
            response = requests.get(url)
            response.raise_for_status()  # Raise an error for bad status codes
            soup = BeautifulSoup(response.content, 'html.parser')

            # Extract description
            description_row = soup.find('td', text="Description:")
            description = description_row.find_next_sibling('td').text.strip() if description_row else None

            # Extract category
            category_row = soup.find('td', text="Category:")
            category = []
            if (category_row):
                category = [a.text for a in category_row.find_next_sibling('td').find_all('a')]

            # Extract value type
            value_type_row = soup.find('td', text="Value Type")
            value_type = value_type_row.find_next_sibling('td').text.strip() if value_type_row else None

            # Extract item count
            item_count_row = soup.find('a', href="help.cgi?cd=item_count")
            item_count = item_count_row.find_next('td', class_='int_blu').text.strip() if item_count_row else None

            # Extract item type
            item_type_row = soup.find('a', href="help.cgi?cd=item_type")
            item_type = item_type_row.find_next('td', class_='txt_blu').text.strip() if item_type_row else None

            # Extract stability
            stability_row = soup.find('a', href="help.cgi?cd=stability")
            stability = stability_row.find_next('td', class_='txt_blu').text.strip() if stability_row else None

            # Extract strata
            strata_row = soup.find('a', href="help.cgi?cd=strata")
            strata = strata_row.find_next('td', class_='txt_blu').text.strip() if strata_row else None

            # Extract array
            array_row = soup.find('a', href="help.cgi?cd=array")
            array = array_row.find_next('td', class_='txt_blu').text.strip() if array_row else None
            
            # Extract number of participants
            participants_row = soup.find('a', href="help.cgi?cd=participant")
            participants = participants_row.find_next('td', class_='int_blu').text.strip() if participants_row else None

            metadata = {
                "Description": description,
                "Category": category,
                "Value Type": value_type,
                "Item Count": item_count,
                "Item Type": item_type,
                "Stability": stability,
                "Strata": strata,
                "Array": array,
                "Participants": participants,
                "Url": url,
            }

            return metadata

        except requests.RequestException as e:
            print(f"Attempt {attempt + 1} failed: {e}")
            time.sleep(time_duration)
    
    return None

def get_phenotype_data_for_field(phenotype_df, field_id, verbose=False):
    # Convert field_id to string and append the specific delimiter, e.g., a period
    prefix = str(field_id) + "."
    if verbose:
        print("Prefix:", prefix)
    # Filter columns based on this prefix
    columns = [col for col in phenotype_df.columns if col.startswith(prefix)]
    if verbose:
        print("Filtered columns:", columns)
    return phenotype_df[columns]

def get_phenotype_data_for_field(phenotype_df, field_id, verbose=False, id_column=None):
    # Convert field_id to string and append the specific delimiter, e.g., a period
    prefix = str(field_id) + "."
    if verbose:
        print("Prefix:", prefix)
    
    # Filter columns based on this prefix
    columns = [col for col in phenotype_df.columns if col.startswith(prefix)]
    if verbose:
        print("Filtered columns:", columns)
    
    # If id_column is provided and exists in the DataFrame, include it in the result
    if id_column and id_column in phenotype_df.columns:
        columns = [id_column] + columns
    
    # Select the filtered columns
    filtered_df = phenotype_df[columns]
    
    # Set the id_column as the index if provided
    if id_column and id_column in filtered_df.columns:
        filtered_df = filtered_df.set_index(id_column)
    
    return filtered_df

def extract_field_data(ukb_phenotypes_df, field_ids):
    df = pd.DataFrame()
    executed_fields = []

    for field_id in field_ids:
        tmp = get_phenotype_data_for_field(ukb_phenotypes_df, field_id=field_id, verbose=False)
        print(f"Field ID: {field_id}")
        # print(f"Shape: {tmp.shape}")
        print("Head:")
        print(tmp.head())

        # calculate n columns
        n_columns = tmp.shape[1]
        first_column = tmp.columns[0]
        try:
            # n_instances = len(tmp.columns.str.split('.').str[1].unique())
            n_instances = len(pd.Index(tmp.columns).str.split('.').str[1].unique())
        except Exception:
            n_instances = 0

        try:
            n_cases =  len(pd.Index(tmp.columns).str.split('.').str[2].unique())
        except Exception:
            n_cases = 0

        unique_values = pd.unique(tmp.to_pandas().values.ravel())
        unique_set = set(unique_values)  # Convert list of unique values to a set for efficient membership testing
        print(f"Unique values: {unique_values}")
        metadata = get_biobank_metadata(field_id, number_attempts=10, time_duration=5)
        print(f"Metadata: {metadata}")

        data_type = tmp[first_column].dtype
        # print(f"Data type: {data_type}")
        
        contains_binary = set(pd.Series(unique_values).dropna().tolist()) == {-3, -1, 0, 1}
        
        # assesment that unique_values contains only 1, 0, -1, and -3 values, avoid NaN values
        # contains_binary = {-3, -1, 0, 1} == set(unique_set)  # Check if all -3, -1, 0, and 1 are exactly in the set
        print(f"Contains binary: {contains_binary}")
        
        contains_yes_no = {'Yes', 'No'}.issubset(unique_set)  # Check if both 'Yes' and 'No' are exactly in the set
        
        # check that values it is possible convert data to datetime format, avoid NaN values, and add the column date_format to the dataframe
        if data_type == 'object':
            try:
                tmp[first_column] = pd.to_datetime(tmp[first_column])
                # print("Data converted to datetime")
                data_format = True
            except:
                # print("Data not converted to datetime")
                data_format = False
        else:
            # print("Data not converted to datetime")
            data_format = False
        
        #######
        if data_type == 'object' and contains_yes_no:
            # print("Yes/No values found")
            df = df._append({
                'field_id': field_id, 
                'data_type': str(data_type),
                'number_uniq_values': len(unique_values),
                'number_columns': n_columns,
                'number_instances': n_instances,
                'number_cases': n_cases,
                'binary': contains_binary,
                # 'Yes_No': contains_yes_no,
                'date_format': data_format,
                'description': metadata['Description'],
                'category': metadata['Category'],
                'value_type': metadata['Value Type'],
                'participants': metadata['Participants'],
                'item_count': metadata['Item Count'],
                'item_type': metadata['Item Type'],
                'stability': metadata['Stability'],
                'strata': metadata['Strata'],
                'array': metadata['Array'],
                'url': metadata['Url'],
            }, ignore_index=True)
        else:
            # print("Yes/No values not found")
            df = df._append({
                'field_id': field_id, 
                'data_type': str(data_type),
                'number_uniq_values': len(unique_values),
                'number_columns': n_columns,
                'number_instances': n_instances,
                'number_cases': n_cases,
                # 'Yes_No': contains_yes_no,
                'contains_binary': contains_binary,
                'date_format': data_format,
                'description': metadata['Description'],
                'category': metadata['Category'],
                'value_type': metadata['Value Type'],
                'participants': metadata['Participants'],
                'item_count': metadata['Item Count'],
                'item_type': metadata['Item Type'],
                'stability': metadata['Stability'],
                'strata': metadata['Strata'],
                'array': metadata['Array'],
                'url': metadata['Url'],
            }, ignore_index=True)

        executed_fields.append(field_id)
        total_fields = len(field_ids)
        executed_count = len(executed_fields)
        print(f"Executed {executed_count}/{total_fields} fields")
        
        print("")

    return df

################################################################################
# function to transform data
################################################################################
def calculate_mean_mean_df(ukb_phenotypes_df, field_id):
    # Get data for the specific field_id
    data_field_id = get_phenotype_data_for_field(ukb_phenotypes_df, field_id=field_id)
    
    # Determine the number of instances
    instances = [element.split('.')[1] for element in data_field_id.columns]

    # Create an empty dataframe
    df = pd.DataFrame()

    # Process each instance
    for instance in instances:
        prefix = f"{field_id}.{instance}."
        
        # Get columns with the prefix
        columns = [col for col in data_field_id.columns if col.startswith(prefix)]

        data_instance = data_field_id[columns]

        # Calculate mean value by row
        mean_value_by_row = data_instance.to_pandas().mean(axis=1)
        
        # Append the mean values to the dataframe
        df = df._append(mean_value_by_row, ignore_index=True)
    
    # Transpose the dataframe
    df = df.T
    
    # Calculate mean of means
    mean_mean_df = df.mean(axis=1)
    mean_mean_df = pd.DataFrame(mean_mean_df)

    # Set the column name using field_id
    mean_mean_df.columns = [field_id]

    # Convert to Polars DataFrame
    mean_mean_df = pl.DataFrame(mean_mean_df)
    
    return mean_mean_df



def calculate_mean_mean_df(ukb_phenotypes_df, field_id):
    # Get data for the specific field_id
    data_field_id = get_phenotype_data_for_field(ukb_phenotypes_df, field_id=field_id)

    # Determine the unique instance numbers
    instances = set(element.split('.')[1] for element in data_field_id.columns)
    
    # Initialize a list to collect mean values for each instance
    mean_values = []

    # Process each instance
    for idx, instance in enumerate(instances):
        prefix = f"{field_id}.{instance}."
        
        # Select columns with the prefix
        columns = [col for col in data_field_id.columns if col.startswith(prefix)]

        # Select data for the current instance
        data_instance = data_field_id.select(columns)
        
        # Calculate mean value by row using Polars expressions
        mean_value_by_row = data_instance.with_columns(pl.mean_horizontal(pl.all()).alias(f'row_mean_{idx}'))

        # Collect the mean values
        mean_values.append(mean_value_by_row[f'row_mean_{idx}'])

    # Combine all mean values into a single DataFrame
    df_means = pl.DataFrame(mean_values)

    # # Calculate the mean of means across rows
    mean_mean_df = df_means.with_columns(pl.mean_horizontal(pl.all()).alias(f'{field_id}'))
    
    # select the column with field_id name
    mean_mean_df = pl.DataFrame(mean_mean_df[f'{field_id}']) 
    
    return mean_mean_df

def calculate_means_for_multiple_fields(ukb_phenotypes_df, field_ids):
    result_dfs = []

    # Iterate over each field_id and calculate the mean_mean_df
    for field_id in field_ids:
        mean_mean_df = calculate_mean_mean_df(ukb_phenotypes_df, field_id)
        result_dfs.append(mean_mean_df)

    # Concatenate all result DataFrames horizontally
    final_df = pl.concat(result_dfs, how='horizontal')

    return final_df


################################################################################
def calculate_max_one_hot_data(ukb_phenotypes_df, field_id):
    # Retrieve field data
    field_data = get_phenotype_data_for_field(ukb_phenotypes_df, field_id)
    
    # One-hot encode the field data
    field_data_one_hot = field_data.to_dummies()
    
    # Remove columns that end with 'null'
    columns_to_remove = [col for col in field_data_one_hot.columns if col.endswith('null')]
    field_data_one_hot = field_data_one_hot.drop(columns_to_remove)
    
    # Extract unique coding values from the column names
    coding = pd.Index(field_data_one_hot.columns).str.split('_').str[1].unique().tolist()
    
    # Remove negative values from the coding list
    # coding = [code for code in coding if float(code) >= 0]
    
    print(coding)
    
    new_coding = []
    
    for code in coding:
        try:
            if float(code) > 0:
                new_coding.append(code)
        except ValueError:
            # If conversion is not possible, classify as good
            new_coding.append(code)
        
    coding = new_coding
    
    max_values = []
    
    for code in coding:
        code_columns = [col for col in field_data_one_hot.columns if col.endswith(code)]
        data_code = field_data_one_hot[code_columns]
        
        max_data_code = data_code.with_columns(pl.max_horizontal(pl.all()).alias(f'{field_id}_{code}_max'))
        
        max_values.append(max_data_code[f'{field_id}_{code}_max'])
    
    # Combine the max values into a DataFrame
    df_max = pl.DataFrame(max_values)
    df_max.columns = [col.replace('_max', '') for col in df_max.columns]
    
    return df_max

def calculate_max_one_hot_for_multiple_fields(ukb_phenotypes_df, field_ids):
    result_dfs = []
    index = 0

    # Iterate over each field_id and calculate the mean_mean_df
    for field_id in field_ids:
        index += 1
        print(f"Processing field {field_id} | Index: {index} of {len(field_ids)}")
        
        max_one_hot_data = calculate_max_one_hot_data(ukb_phenotypes_df, field_id)
        result_dfs.append(max_one_hot_data)

    # Concatenate all result DataFrames horizontally
    final_df = pl.concat(result_dfs, how='horizontal')

    return final_df


################################################################################
os.environ['POLARS_MAX_THREADS'] = '20'

ukb_phenotypes_df = pl.read_csv(
    'raw/ukb-phenotypes/ukb673991.tab',
    separator='\t',
    infer_schema_length=1000000,  # Increase schema inference length
    # dtypes={'f.630.0.0': pl.Int64},  # Specify data type for specific columns if known
    ignore_errors=False,  # Ignore parsing errors
    null_values='NA'  # Treat 'NA' as null values
)

# remove f. in the column names of ukb_phenotypes_df
ukb_phenotypes_df.columns = [col.replace('f.', '') for col in ukb_phenotypes_df.columns]

# read the RealDiseaseCodesForPatients.txt file
# ukb_phenotypes_df = pd.read_csv('data/UKB-phenotypes/ukb-phenotypes-merged.tsv', sep='\t')
ukb_phenotypes_df.head()
ukb_phenotypes_df.shape

# select column names, and select first part of the column name befor the '.', and get unique values
field_ids = np.unique([col.split(".")[0] for col in ukb_phenotypes_df.columns]).tolist()
field_ids.remove('eid')
field_ids
len(field_ids)


ali_field_analysis = [
    "34", "48", "49", "137", "4079", "4080", "12144", 
    "21001", "21002", "30000", "30010", "30040", "30050", 
    "30080", "30100", "30130", "40008", "1259", "20551", 
    "2492", "6153", "6156", "1558", "2443", "4041", "6154", 
    "6177", "1618", "2453", "6155", "41270"
]

################################################################################
# assesment of transformation of phenotype data
################################################################################

df2 = extract_field_data(ukb_phenotypes_df, field_ids)

# save df to tsv file, with columns names, without index
df2.to_csv('results/google-drive/preparing-phenotypes/ukb-phenotypes-processed-raw2.tsv', sep='\t', index=False)

# read the ukb-phenotypes-processed-raw.tsv file
df = pd.read_csv('results/google-drive/preparing-phenotypes/ukb-phenotypes-processed-raw2.tsv', sep='\t')
df['value_type_simple'] = df['value_type'].str.split(',').str[0]
df['array_simple'] = df['array'].str.split(' ').str[0]

# remove coma from item_count column
df['item_count'] = df['item_count'].str.replace(',', '')
df.item_count = df.item_count.astype(int)

# remove from df value_type equal "Date" or "Time"
df = df[~df['value_type'].str.contains('Date')]
df = df[~df['value_type'].str.contains('Time')]
df = df[~df['value_type'].str.contains('Text')]
df = df[~df['value_type'].str.contains('Continuous, years')]
df = df[~df['value_type'].str.contains('Compound')]
df = df[df['number_instances'] <= 5]

# replace ',' with  '' in participants column
df['participants'] = df['participants'].str.replace(',', '')
df.participants = df.participants.astype(int)
df[df['participants'] < 5000]


# create new column name assessment_data, at default value is manual
df['assessment_data'] = 'manual'
# for Integer values in values_type_simple, set assessment_data to mean_mean
df.loc[df['value_type_simple'] == 'Integer', 'assessment_data'] = 'mean_mean'

# for Continuous values in values_type_simple, set assessment_data to mean_mean
df.loc[df['value_type_simple'] == 'Continuous', 'assessment_data'] = 'mean_mean'

# for Categorical (single) values in values_type_simple, set assessment_data to one_hot_max
df.loc[df['value_type_simple'] == 'Categorical (single)', 'assessment_data'] = 'one_hot_max'

# for Categorical (multiple) values in values_type_simple, set assessment_data to one_hot_max
df.loc[df['value_type_simple'] == 'Categorical (multiple)', 'assessment_data'] = 'one_hot_max'

# for containing binary values in contains_binary, set assessment_data to binary_max
df.loc[df['contains_binary'] == True, 'assessment_data'] = 'binary_max'

# save df to tsv file, with columns names, without index
df.to_csv('results/google-drive/preparing-phenotypes/ukb-phenotypes-assessed-transformation.tsv', sep='\t', index=False)

################################################################################
# prepare function to calculate the max value by each row in the dataframe
################################################################################

# filter mean_mean in assessment_data column
df_mean_mean = df[df['assessment_data'] == 'mean_mean']

df_processing_mean_mean = calculate_means_for_multiple_fields(ukb_phenotypes_df, df_mean_mean['field_id'].tolist())


################################################################################
# prepare functionj one hot encoding
################################################################################



# get_phenotype_data_for_field(ukb_phenotypes_df, field_id=54, verbose=True, id_column='eid')  

# # minus values convert to null values
# tmp = get_phenotype_data_for_field(ukb_phenotypes_df, field_id=54)
# tmp.head()

# tmp_one_hot = tmp.to_dummies()

# # filter elemnt of list which ended with 'null'
# columns_to_remove = [col for col in tmp_one_hot.columns if col.endswith('null')]

# tmp_one_hot = tmp_one_hot.drop(columns_to_remove)

# tmp_one_hot.head()



# field_id = 54

# field_data = get_phenotype_data_for_field(ukb_phenotypes_df, field_id=54)

# field_data_one_hot = field_data.to_dummies()

# columns_to_remove = [col for col in field_data_one_hot.columns if col.endswith('null')]

# field_data_one_hot = field_data_one_hot.drop(columns_to_remove)

# # select second part of the column name after the '_', and get unique values
# coding = pd.Index(field_data_one_hot.columns).str.split('_').str[1].unique().tolist()

# # add to list '-1' element
# coding.append('-1')

# coding.append('0')
# coding.append('0.5')
# coding.append('<1')
# # from list remove lower than 0 values
# coding = [code for code in coding if float(code) >= 0]
# coding
# coding
# max_values = []

# for code in coding:
#     code_columns = [col for col in field_data_one_hot.columns if col.endswith(code)]
    
#     print(code_columns)
    
#     data_code = field_data_one_hot[code_columns].head()
    
#     max_data_code = data_code.with_columns(pl.max_horizontal(pl.all()).alias(f'{field_id}_{code}_max'))
    
#     print(max_data_code.head())
    
#     max_values.append(max_data_code[f'{field_id}_{code}_max'])
    
    
# print(max_values)

# df_max = pl.DataFrame(max_values)
# df_max.columns = [col.replace('_max', '') for col in df_max.columns]
# df_max


# filter in assessment_data column using list of values
df_one_hot_max = df[df['assessment_data'].isin(['one_hot_max', 'binary_max'])]


# from df_one_hot_max remove filed_id for '400002', '41202'
icd10_field_id_to_remove = [40002, 41202, 41204, 41201, 40006, 40001]
icd9_field_it_to_remove = [40008, 40012, 40019, 40021, 40005, 40011, 40004, 40006, 40013]


df_one_hot_max = df_one_hot_max[~df_one_hot_max['field_id'].isin(icd10_field_id_to_remove)]
df_one_hot_max = df_one_hot_max[~df_one_hot_max['field_id'].isin(icd9_field_it_to_remove)]


calculate_max_one_hot_data(ukb_phenotypes_df, field_id=54)

df_processing_one_hot_max = calculate_max_one_hot_for_multiple_fields(ukb_phenotypes_df, field_ids = df_one_hot_max['field_id'].tolist())

df_processing_one_hot_max.head()

# combine df_mean_mean and df_processing_one_hot_max, 
df_processing = pl.concat([df_processing_mean_mean, df_processing_one_hot_max], how='horizontal')

df_processing

# add eid column to the df_processing, df_processing has no eid column
df_processing = pl.concat([pl.DataFrame(ukb_phenotypes_df['eid']), df_processing], how='horizontal')

df_processing_mean_mean.shape
df_processing_one_hot_max.shape
df_processing.shape

# save df_processing to tsv file, with columns names, without index
df_processing.write_csv('data/UKB-phenotypes/ukb-phenotypes-processed-raw.tsv', separator='\t')

del ukb_phenotypes_df

# perform imputation of missing values in the df_processing
df_processing_imputed = impute_missing_values(df_processing)




################################################################################
# read the ukb-phenotypes-processed-raw.tsv file using polars
################################################################################
os.environ['POLARS_MAX_THREADS'] = '20'

ukb_phenotypes_df = pl.read_csv(
    'data/UKB-phenotypes/ukb-phenotypes-processed-raw.tsv',
    separator=',',
    infer_schema_length=1000000,  # Increase schema inference length
    # dtypes={'f.630.0.0': pl.Int64},  # Specify data type for specific columns if known
    ignore_errors=False,  # Ignore parsing errors
    null_values='NA'  # Treat 'NA' as null values
)



from tqdm import tqdm

def standardize_columns(df):
    """
    Standardize each column of the Polars DataFrame separately.

    Parameters:
    df (pl.DataFrame): The DataFrame to standardize.

    Returns:
    pl.DataFrame: The standardized DataFrame.
    """
    standardized_data = {}
    
    # Use tqdm to show a progress bar
    for column in tqdm(df.columns, desc="Standardizing columns"):
        mean = df[column].mean()
        std_dev = df[column].std()
        standardized_data[column] = (df[column] - mean) / std_dev
    
    standardized_df = pl.DataFrame(standardized_data)
    
    return standardized_df

def impute_missing_values(df: pl.DataFrame) -> pl.DataFrame:
    """
    Impute missing values in each column separately with the column mean.
    
    Parameters:
    df (pl.DataFrame): The input Polars DataFrame with missing values.
    
    Returns:
    pl.DataFrame: A new Polars DataFrame with missing values imputed.
    """
    # Iterate over each column
    for column in df.columns:
        # Calculate the mean of the column excluding missing values
        mean_value = df[column].mean()
        
        # Replace missing values with the mean
        df = df.with_columns(
            pl.when(pl.col(column).is_null())
            .then(mean_value)
            .otherwise(pl.col(column))
            .alias(column)
        )
    
    return df
def impute_missing_values(df: pl.DataFrame) -> pl.DataFrame:
    """
    Impute missing values (None and NaN) in each column separately with the column mean.
    
    Parameters:
    df (pl.DataFrame): The input Polars DataFrame with missing values.
    
    Returns:
    pl.DataFrame: A new Polars DataFrame with missing values imputed.
    """
    # Iterate over each column
    for column in df.columns:
        # Calculate the mean of the column excluding missing values
        mean_value = df.filter(pl.col(column).is_not_null() & pl.col(column).is_finite())[column].mean()
        
        # Replace None and NaN values with the mean
        df = df.with_columns(
            pl.when(pl.col(column).is_null() | pl.col(column).is_nan())
            .then(mean_value)
            .otherwise(pl.col(column))
            .alias(column)
        )
    
    return df

df_raw_standardized = standardize_columns(df_processing[:, 1:10000])
df_raw_standardized = impute_missing_values(df_raw_standardized)

# check that df_raw_standardized has NaN values





# performing pca on the df_raw_standardized
from sklearn.decomposition import PCA
import numpy as np

import matplotlib.pyplot as plt


data_np = df_raw_standardized.to_numpy()
pca = PCA(n_components=2)
principalComponents = pca.fit_transform(data_np)
