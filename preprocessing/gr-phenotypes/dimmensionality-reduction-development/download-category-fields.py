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
import requests
from bs4 import BeautifulSoup
import multiprocessing


def get_biobank_metadata(field_id, number_attempts=5, time_duration=5):
    """
    Function to download the metadata information from the UK Biobank website with retries.

    Parameters:
    field_id (int): The field ID to download the metadata.
    number_attempts (int): The number of attempts to make to download the metadata.
    time_duration (int): The time duration to wait between attempts in seconds.

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

            # Determine last category
            last_category = category[-1] if category else None

            metadata = {
                "field_id": field_id,
                "description": description,
                "category": category,
                "last_category": last_category
            }

            return metadata

        except requests.RequestException as e:
            print(f"Attempt {attempt + 1} failed: {e}")
            time.sleep(time_duration)
    
    return None

def worker(field_id):
    return get_biobank_metadata(field_id)

# read the field_ids from the file
field_ids = pd.read_csv('results/google-drive/preparing-phenotypes/field.tsv', sep='\t') 

# fielter the fields id below 13000
field_ids = field_ids[field_ids['field_id'] < 130000]
field_id = field_ids.field_id.tolist()
field_id = [str(i) for i in field_id]


num_processes = 10
# Create a pool of workers
with multiprocessing.Pool(processes=num_processes) as pool:
    results = list(tqdm(pool.imap(worker, field_id), total=len(field_id)))



# convert dictionary to dataframe with three columns: field_id, categories, last_category, categories as a list
field_ids_categories = pd.DataFrame(results)
field_ids_categories['categories'] = field_ids_categories['category'].apply(lambda x: ', '.join(x))
field_ids_categories = field_ids_categories.drop(columns=['category'])
field_ids_categories = field_ids_categories[['field_id', 'description', 'categories', 'last_category']]

field_ids_categories.head()

field_ids_categories.to_csv('results/google-drive/preparing-phenotypes/fields_categories.tsv', sep='\t', index=False)