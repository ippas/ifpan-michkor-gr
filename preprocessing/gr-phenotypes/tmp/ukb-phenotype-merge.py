# import modules
import os
import pandas as pd

print(os.getcwd())

# find tsv files in the directory
tsv_files = [file for file in os.listdir('raw/dataset/') if file.endswith('.tsv')]

print(tsv_files)

# read the first file
df_1 = pd.read_csv('raw/dataset/ukb51759.tsv', sep='\t')
df_2 = pd.read_csv('raw/dataset/ukb49742.tsv', sep='\t')
df_3 = pd.read_csv('raw/dataset/ukb48815.tsv', sep='\t')

# check dimensions of the dataframe
print(df_1.shape)
print(df_2.shape)
print(df_3.shape)

# select column names for each dataframe
print(df_1.columns)
print(df_2.columns)
print(df_3.columns)

# select f.eid and convert to list
eid_1 = df_1[['f.eid']].iloc[:, 0].tolist()
eid_2 = df_2[['f.eid']].iloc[:, 0].tolist()
eid_3 = df_3[['f.eid']].iloc[:, 0].tolist()

# get common elements between three lists
common_eid = list(set(eid_1) & set(eid_2) & set(eid_3))

# length of the all four list
print(len(eid_1))
print(len(eid_2))
print(len(eid_3))
print(len(common_eid))

# filter df_1, df_2, df_3 for common_eid
df_1 = df_1[df_1['f.eid'].isin(common_eid)]
df_2 = df_2[df_2['f.eid'].isin(common_eid)]
df_3 = df_3[df_3['f.eid'].isin(common_eid)]


# find common columns between df_2 and df_3 and remove f.eid
common_columns = df_2.columns.intersection(df_3.columns).difference(['f.eid'])
print(common_columns)
print(len(common_columns))

# print df_2 for columns f.eid and 'f.87.0.0'
print(df_2[['f.eid', 'f.87.0.0']].tail())   
print(df_3[['f.eid', 'f.87.0.0']].tail())   

# from df_3 remove common_columns
df_3 = df_3.drop(columns=common_columns)

# merge df_1, df_2, df_3 to one dataframe
ukb_phenotypes_df = pd.merge(df_1, df_2, on='f.eid', how='outer')
ukb_phenotypes_df = pd.merge(ukb_phenotypes_df, df_3, on='f.eid', how='outer')

print(ukb_phenotypes_df.shape)

# remove 'f.' prefix from column names
ukb_phenotypes_df.columns = ukb_phenotypes_df.columns.str.replace('f.', '')

# save the dataframe to a file
ukb_phenotypes_df.to_csv('data/UKB-phenotypes/ukb-phenotypes-merged.tsv', sep='\t', index=False)

########################################################################################
# Read the output_ukb.txt file
########################################################################################

# Read only the first row to get the column names
output_ukb = pd.read_csv('/net/pr2/projects/plgrid/plggneuromol/Alireza/UkbPhenotypes/output_ukb.txt', sep='\t', nrows=0)

# in each column name for output_ukb, replace '-' with '_'
output_ukb.columns = output_ukb.columns.str.replace('-', '.')

# The column names are now stored in output_ukb.columns
print(output_ukb.columns)

# find the common columns between ukb_phenotypes_df and output_ukb
print(ukb_phenotypes_df.columns.intersection(output_ukb.columns))


# find columns which occur only in ukb_phenotypes_df
only_in_ukb_phenotypes = ukb_phenotypes_df.columns.difference(output_ukb.columns)
print(only_in_ukb_phenotypes)

# find columns which occur only in output_ukb
only_in_output_ukb = output_ukb.columns.difference(ukb_phenotypes_df.columns)
print(only_in_output_ukb.to_list())


