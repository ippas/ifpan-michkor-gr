# import modules
import os
import pandas as pd
import numpy as np
import requests
from bs4 import BeautifulSoup
import time
import openpyxl
import multiprocessing as mp
import polars as pl
from concurrent.futures import ThreadPoolExecutor, as_completed
import concurrent.futures
from tqdm import tqdm


################################################################################
# fucntions
################################################################################

################################################################################
# function to download data from UK Biobank
################################################################################

# function to download the UK Biobank phenotypes
def get_biobank_field_showcase(field_id, number_attempts=5, time_duration=5):
    """
    Function to download the metadata information from the UK Biobank website.
    
    Parameters:
    field_id (int): The field ID to download the metadata.
    number_attempts (int): The number of attempts to make to download the metadata.
    time_duration (int): The time duration to wait between attempts.
    
    Returns:
    dict: The metadata information in a dictionary.
    """
    
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

            # Extract sexed information
            sexed_row = soup.find('a', href="help.cgi?cd=sexed")
            sexed = sexed_row.find_next('td', class_='txt_blu').text.strip() if sexed_row else None

            # Extract instances information
            instances_row = soup.find('a', href="help.cgi?cd=instances")
            instances = instances_row.find_next('td', class_='txt_blu').text.strip() if instances_row else None

            # Extract data coding number
            data_section = soup.find('div', class_='tabbertab').find('h2', text='Data')
            data_coding = None
            
            # Extract debut information
            debut_row = soup.find('a', href="help.cgi?cd=debut")
            debut = debut_row.find_next('td', class_='txt_blu').text.strip() if debut_row else None
            
            # Extract version information
            version_row = soup.find('a', href="help.cgi?cd=version")
            version = version_row.find_next('td', class_='txt_blu').text.strip() if version_row else None
            
            # Extract cost tier information
            cost_tier_row = soup.find('a', href="help.cgi?cd=cost_tier")
            cost_tier = cost_tier_row.find_next('td', class_='txt_blu').text.strip() if cost_tier_row else None

            if data_section:
                data_text = data_section.find_next_sibling(text=True)
                if data_text and 'Data-Coding' in data_text:
                    data_coding_tag = data_section.find_next('a', class_='basic', href=True)
                    if data_coding_tag:
                        href = data_coding_tag['href']
                        data_coding = href.split('=')[-1]

            metadata = {
                "description": description,
                "category": category,
                "participants": participants,
                "item_count": item_count,
                "value_type": value_type,
                "stability": stability,
                "value_type": value_type, # is not prepared
                "item_type": item_type,
                "strata": strata,
                "sexed": sexed,
                "instances": instances,
                "array": array,
                "debut": debut, # is not prepared
                "version": version, # is not prepared
                "cost_tier": cost_tier, # is not prepared
                "data_coding": data_coding,
                "url": url,
            }

            return metadata

        except requests.RequestException as e:
            print(f"Attempt {attempt + 1} failed: {e}")
            time.sleep(time_duration)
    
    return None

# function to get meaning of data coding
def get_data_coding_ukb(data_coding, number_attempts=5, time_duration=5):
    # description of function
    """
    Function to download the data coding information from the UK Biobank website.
    
    Parameters:
    data_coding (str): The data coding value to download.
    number_attempts (int): The number of attempts to make to download the data coding.
    time_duration (int): The time duration to wait between attempts.
    
    Returns:
    pd.DataFrame: The data coding information in a DataFrame.
    """
    
    # Define specific messages for certain data_coding values
    specific_messages = {
        '19': 'There are ICD10 coding.',
        '2': 'There are SOC2000 coding.',
        '87': 'There are ICD9 coding.',
    }

    # Define different URL structure for certain data_coding values
    alternative_url_codings = {'123', '240', '259', '4', '744'}

    # Check if the data_coding value has a specific message
    if str(data_coding) in specific_messages:
        print(specific_messages[str(data_coding)])
        return None
    
    # Construct the URL with the given data_coding
    if str(data_coding) in alternative_url_codings:
        url = f"https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id={data_coding}&nl=1"
    else:
        url = f"https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id={data_coding}"
    
    for attempt in range(number_attempts):
        try:
            response = requests.get(url)
            response.raise_for_status()  # Raise an error for bad status codes
            page_content = response.content

            # Parse the HTML
            soup = BeautifulSoup(page_content, 'html.parser')

            # Find the table
            table = soup.find('table', {'border': '0', 'cellspacing': '2'})

            # Extract table headers
            headers = [header.text for header in table.find_all('th')]

            # Extract table rows
            rows = []
            for row in table.find_all('tr')[1:]:  # Skip header row
                cols = row.find_all('td')
                cols = [ele.text.strip() for ele in cols]
                rows.append(cols)

            # Create a DataFrame
            df = pd.DataFrame(rows, columns=headers)
            
            return df

        except requests.RequestException as e:
            print(f"Attempt {attempt + 1} failed: {e}")
            time.sleep(time_duration)
    
    return None

# function to get data coding for multiple data coding values
def get_biobank_multiple_data_coding(data_coding_ukb):
    """
    Downloads data coding for all the data_coding values provided.

    Parameters:
    data_coding_ukb (list): List of data coding identifiers.

    Returns:
    pd.DataFrame: DataFrame containing all the combined data coding information.
    """
    data_coding_ukb_all = {}
    number = 0

    for data_coding in data_coding_ukb:
        number += 1
        print(f"Downloading data coding {number}/{len(data_coding_ukb)}")
        print(data_coding)
        
        data_coding_ukb_all[data_coding] = get_data_coding_ukb(data_coding)
        
    print("Combining data coding")
     
    data_coding_ukb_all_df = pd.concat(data_coding_ukb_all).reset_index(level=1, drop=True).reset_index()
    data_coding_ukb_all_df = data_coding_ukb_all_df.iloc[:, 0:3]
    data_coding_ukb_all_df.columns = ['data_coding', 'coding', 'meaning']
    
    return data_coding_ukb_all_df


################################################################################
# function to processing local data from UK Biobank
################################################################################

# function to load the UK Biobank phenotypes
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

# def extract_field_data(ukb_phenotypes_df, field_ids):
#     df = pd.DataFrame()
#     executed_fields = []

#     for field_id in field_ids:
#         tmp = get_phenotype_data_for_field(ukb_phenotypes_df, field_id=field_id, verbose=False)
#         print(f"Field ID: {field_id}")
#         # print(f"Shape: {tmp.shape}")
#         print("Head:")
#         print(tmp.head())

#         # calculate n columns
#         n_columns = tmp.shape[1]
#         first_column = tmp.columns[0]
#         try:
#             # n_instances = len(tmp.columns.str.split('.').str[1].unique())
#             n_instances = len(pd.Index(tmp.columns).str.split('.').str[1].unique())
#         except Exception:
#             n_instances = 0

#         try:
#             n_cases =  len(pd.Index(tmp.columns).str.split('.').str[2].unique())
#         except Exception:
#             n_cases = 0

#         unique_values = pd.unique(tmp.to_pandas().values.ravel())
#         unique_set = set(unique_values)  # Convert list of unique values to a set for efficient membership testing
#         print(f"Unique values: {unique_values}")
#         metadata = get_biobank_field_showcase(field_id, number_attempts=10, time_duration=5)
#         print(f"Metadata: {metadata}")

#         data_type = tmp[first_column].dtype
#         # print(f"Data type: {data_type}")
        
#         contains_binary = set(pd.Series(unique_values).dropna().tolist()) == {-3, -1, 0, 1}
        
#         # assesment that unique_values contains only 1, 0, -1, and -3 values, avoid NaN values
#         # contains_binary = {-3, -1, 0, 1} == set(unique_set)  # Check if all -3, -1, 0, and 1 are exactly in the set
#         print(f"Contains binary: {contains_binary}")
        
#         contains_yes_no = {'Yes', 'No'}.issubset(unique_set)  # Check if both 'Yes' and 'No' are exactly in the set
        
#         # check that values it is possible convert data to datetime format, avoid NaN values, and add the column date_format to the dataframe
#         if data_type == 'object':
#             try:
#                 tmp[first_column] = pd.to_datetime(tmp[first_column])
#                 # print("Data converted to datetime")
#                 data_format = True
#             except:
#                 # print("Data not converted to datetime")
#                 data_format = False
#         else:
#             # print("Data not converted to datetime")
#             data_format = False
        
#         #######
#         if data_type == 'object' and contains_yes_no:
#             # print("Yes/No values found")
#             df = df._append({
#                 'field_id': field_id, 
#                 'data_type': str(data_type),
#                 'number_uniq_values': len(unique_values),
#                 'number_columns': n_columns,
#                 'number_instances': n_instances,
#                 'number_cases': n_cases,
#                 'binary': contains_binary,
#                 # 'Yes_No': contains_yes_no,
#                 'date_format': data_format,
#                 'description': metadata['Description'],
#                 'category': metadata['Category'],
#                 'value_type': metadata['Value Type'],
#                 'participants': metadata['Participants'],
#                 'item_count': metadata['Item Count'],
#                 'item_type': metadata['Item Type'],
#                 'stability': metadata['Stability'],
#                 'strata': metadata['Strata'],
#                 'array': metadata['Array'],
#                 'url': metadata['Url'],
#             }, ignore_index=True)
#         else:
#             # print("Yes/No values not found")
#             df = df._append({
#                 'field_id': field_id, 
#                 'data_type': str(data_type),
#                 'number_uniq_values': len(unique_values),
#                 'number_columns': n_columns,
#                 'number_instances': n_instances,
#                 'number_cases': n_cases,
#                 # 'Yes_No': contains_yes_no,
#                 'contains_binary': contains_binary,
#                 'date_format': data_format,
#                 'description': metadata['Description'],
#                 'category': metadata['Category'],
#                 'value_type': metadata['Value Type'],
#                 'participants': metadata['Participants'],
#                 'item_count': metadata['Item Count'],
#                 'item_type': metadata['Item Type'],
#                 'stability': metadata['Stability'],
#                 'strata': metadata['Strata'],
#                 'array': metadata['Array'],
#                 'url': metadata['Url'],
#             }, ignore_index=True)

#         executed_fields.append(field_id)
#         total_fields = len(field_ids)
#         executed_count = len(executed_fields)
#         print(f"Executed {executed_count}/{total_fields} fields")
        
#         print("")

#     return df

################################################################################
# function to transform data
################################################################################


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
    
    # new_coding = []
    
    # for code in coding:
    #     try:
    #         if float(code) > 0:
    #             new_coding.append(code)
    #     except ValueError:
    #         # If conversion is not possible, classify as good
    #         new_coding.append(code)
        
    # coding = new_coding
    
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


def calculate_mean_max_icd_category(data, field_id):
    """
    Calculate the mean_max value for each unique pattern in the one-hot encoded data based on field_id.
    
    Parameters:
    data: polars.DataFrame
        The data to be processed.
    field_id: str
        The field_id to be processed.
        
    Returns:
    polars.DataFrame
        The mean_max value for each unique pattern.
    """

    # get lengh of pattern
    pattern_length = len(field_id) + 4
    
    # Generate the one-hot encoded data based on field_id
    df_one_hot_max = calculate_max_one_hot_data(ukb_phenotypes_df=data, field_id=field_id)

    # Select the first n characters of the column names
    pattern_columns = [col[:pattern_length] for col in df_one_hot_max.columns]

    # Get unique values of pattern_columns
    unique_patterns = list(set(pattern_columns))

    # Create a dictionary with key as pattern_columns and value as columns with the same pattern
    dict_pattern_columns = {pattern: [col for col in df_one_hot_max.columns if col.startswith(pattern)] for pattern in unique_patterns}

    # Initialize a list to store the pattern and its mean_max value
    max_values = []
    num_patterns = 0

    # Iterate through each pattern and calculate mean_max
    for pattern, columns in dict_pattern_columns.items():
        print(f"Pattern: {pattern}")
        # Select columns with the current pattern
        part_data_pattern = df_one_hot_max.select(columns)
        
        part_data_pattern = pl.DataFrame(part_data_pattern)
        
        # Calculate mean_max for the current pattern
        mean_max_value = part_data_pattern.with_columns(pl.max_horizontal(pl.all()).alias(f'{pattern}_max'))

        # Append the result as a tuple
        max_values.append(pl.DataFrame(mean_max_value[f"{pattern}_max"]))
        
        num_patterns += 1
        # print how many patterns have been processed
        print(f"Processed pattern: {num_patterns}/{len(dict_pattern_columns)})")
        
    # combine series of max_values to a polars dataframe
    result_df = pl.concat(max_values, how='horizontal')
    
    # remove '_max' from the column names
    result_df.columns = [col.replace('_max', '') for col in result_df.columns]

    return result_df

# function to remove columns with prefix from column names
def remove_columns_with_prefix(df, prefix):
    return df.select([col for col in df.columns if not col.startswith(prefix)])