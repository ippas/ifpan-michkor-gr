import sys
import os
import pandas as pd
import numpy as np
import requests
from bs4 import BeautifulSoup
import time
import multiprocessing as mp
import polars as pl
import multiprocessing
from tqdm import tqdm
import requests
sys.path.append('preprocessing/gr-phenotypes/dimmensionality-reduction-development')

from preprocessing_phenotypes_functions import *


########################################################################################
# code for all the field_ids download from DnaNexus
########################################################################################
# read the field_ids from the file
field_ids = pd.read_csv('results/google-drive/preparing-phenotypes/field.tsv', sep='\t') 

# fielter the fields id below 13000
field_ids = field_ids[field_ids['field_id'] < 130000]

field_ids = field_ids.field_id.tolist()

# convert the field_id to string
field_ids = [str(i) for i in field_ids]

########################################################################################
# download the showcase for all the field_ids
########################################################################################
# Download the metadata for all the field_ids
def worker(field_id):
    return field_id, get_biobank_field_showcase(field_id)

# example
get_biobank_field_showcase("41270")

# Download the metadata for all the field_ids
num_processes = 10
# Create a pool of workers
with multiprocessing.Pool(processes=num_processes) as pool:
    results = list(tqdm(pool.imap(worker, field_ids), total=len(field_ids)))

# Combine results into the dictionary
metadata_ukb_all_filed = {field_id: metadata for field_id, metadata in results}

# convert to dataframe
metadata_ukb_all_filed_df = pd.DataFrame(metadata_ukb_all_filed).T
metadata_ukb_all_filed_df = metadata_ukb_all_filed_df.drop(columns=['url'])
metadata_ukb_all_filed_df.head()


########################################################################################
# get unique data_coding values
########################################################################################
metadata_ukb_all_filed_df = pd.read_csv('results/google-drive/preparing-phenotypes/field.tsv', sep='\t') 

data_coding_ukb = metadata_ukb_all_filed_df[metadata_ukb_all_filed_df['encoding_id'].notnull()].encoding_id.unique().tolist()
data_coding_ukb = [str(i) for i in data_coding_ukb]
data_coding_ukb.remove("0")



get_biobank_multiple_data_coding(data_coding_ukb = data_coding_ukb[1:10])




# download the data coding for all the data_coding values
data_coding_ukb_all = {}
number = 0



for data_coding in data_coding_ukb[1:10]:
    number += 1
    print(f"Downloading data coding {number}/{len(data_coding_ukb)}")
    print(data_coding)
    
    data_coding_ukb_all[data_coding] = get_data_coding_ukb(data_coding)


data_coding_ukb_all_df = pd.concat(data_coding_ukb_all).reset_index(level=1, drop=True).reset_index()
# data_coding_ukb_all_df = data_coding_ukb_all_df.drop(columns=['Selectable', 'Node', 'Parent'])
data_coding_ukb_all_df = data_coding_ukb_all_df.iloc[:, 0:3]

data_coding_ukb_all_df.columns = ['data_coding', 'coding', 'meaning']



# save data_coding_ukb_all_df to a tsv file and xlsx file
data_coding_ukb_all_df.to_csv('results/google-drive/preparing-phenotypes/dimmensionality-reduction-development/data-coding-ukb.tsv', sep='\t', index=False)
data_coding_ukb_all_df.to_excel('results/google-drive/preparing-phenotypes/dimmensionality-reduction-development/data-coding-ukb.xlsx', index=False)


# get first 3 columns from df 
data_coding_ukb_all_df = data_coding_ukb_all_df.iloc[:, 0:3]