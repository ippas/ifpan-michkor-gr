import sys
import os
import pandas as pd
import numpy as np
import requests
from bs4 import BeautifulSoup
import time
import multiprocessing as mp

import polars as pl
import tqdm

sys.path.append('preprocessing/gr-phenotypes/dimmensionality-reduction-development')

from preprocessing_phenotypes_functions import *

# fields to analysis
ali_field_analysis = [
    "34", "48", "49", "137", "4079", "4080", "12144", 
    "21001", "21002", "30000", "30010", "30040", "30050", 
    "30080", "30100", "30130", "40008", "1259", "20551", 
    "2492", "6153", "6156", "1558", "2443", "4041", "6154", 
    "6177", "1618", "2453", "6155", "41270"
]

# read feld-id-from-FA.txt as list
with open('data/UKB-phenotypes/dimmensionality-reduction-development/field-id-from-FA.txt', 'r') as file:
    field_id_from_FA = file.read().splitlines()

# read feld-id-from-FA.txt as list
with open('data/UKB-phenotypes/dimmensionality-reduction-development/field-phenocode-from-FA.txt', 'r') as file:
    field_phenocode_from_FA = file.read().splitlines()

# read data
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


# read ukb-phenotypes-assessed-transformation.tsv
df = pd.read_csv(
    'results/google-drive/preparing-phenotypes/ukb-phenotypes-assessed-transformation.tsv',
    sep='\t'
)

############################################
# filter data
############################################
# filter field_id from df by field_id_from_FA
# convert charater list to integer list
field_id_from_FA = list(map(int, field_id_from_FA))
filtered_df =  df[df['field_id'].isin(field_id_from_FA)]
min(filtered_df.participants.tolist())

df_mean_mean = filtered_df[filtered_df['assessment_data'] == 'mean_mean']
df_processing_mean_mean = calculate_means_for_multiple_fields(ukb_phenotypes_df, df_mean_mean['field_id'].tolist())

df_one_hot_max = filtered_df[filtered_df['assessment_data'] == 'one_hot_max']
# filter field_id for 41202
df_one_hot_max[df_one_hot_max['field_id'] == 41202]

df_processing_one_hot_max = calculate_max_one_hot_for_multiple_fields(ukb_phenotypes_df, field_ids = df_one_hot_max['field_id'].tolist())
################################################################################################


################################################################################################
# corrected icd10 code
################################################################################################


df_processing_one_hot_max = remove_columns_with_prefix(df_processing_one_hot_max, '41202_')

icd10_41202 = calculate_mean_max_icd_category(ukb_phenotypes_df, field_id="41202")

# combine all the preprocessing data
df_preprocessing_field_fa = pl.concat([pl.DataFrame(ukb_phenotypes_df['eid']), 
                                       df_processing_mean_mean, 
                                       df_processing_one_hot_max, 
                                       icd10_41202], 
                                      how='horizontal')

################################################################################################
# checking data with field_phenocode_from_FA

len(field_phenocode_from_FA)

df_preprocessing_field_fa.columns

# find in df_preprocessing_field_fa columsn with prefix 2000
[col for col in df_preprocessing_field_fa.columns if col.startswith('2000_')]

[col for col in df_preprocessing_field_fa.columns if col in field_phenocode_from_FA]

# show elements from field_phenocode_from_FA that are not in df_preprocessing_field_fa.columns
not_found_phtnoeytpes = [col for col in field_phenocode_from_FA if col not in df_preprocessing_field_fa.columns]

# filter df by 2000 field id
df[df['field_id'] == 20509]

# find columns with prefix 2000_ in ukb_phenotypes_df
[col for col in ukb_phenotypes_df.columns if col.startswith('981')]

# calculate sum for each column with prefix 20509_
df_preprocessing_field_fa[[col for col in df_processing_one_hot_max.columns if col.startswith('20516_')]].sum()

df_preprocessing_field_fa[[col for col in df_processing_one_hot_max.columns if col.startswith('20509_')]]

ukb_phenotypes_df['20516.0']

################################################################################################

df[df.contains_binary == True].field_id.tolist()

# convert as character
df[df.contains_binary == True].field_id.astype(str).tolist()

[col for col in not_found_phtnoeytpes if col not in df[df.contains_binary == True].field_id.astype(str).tolist()]

# print list, each element in new line
print('\n'.join([col for col in not_found_phtnoeytpes if col not in df[df.contains_binary == True].field_id.astype(str).tolist()]))

# e.g 404 field_id is intiger, why not found in processing data
df[df.field_id == "404"]

ukb_phenotypes_df['404']