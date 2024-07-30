# import modules
import os
import pandas as pd
import numpy as np

print(os.getcwd())

ukb_phenotype_field_ID_list = ['12144', '1259', '137', '1558', '1618', '20118', '20466', '20551', '21001', '21002', '2443', '2453', '2492', '30000', '30010', '30040', '30050', '30080', '30100', '30130', '31', '3159', '34', '40008', '4041', '4079', '4080', '48', '49', '54', '6153', '6154', '6155', '6156', '6177']
print(len(ukb_phenotype_field_ID_list))
print(ukb_phenotype_field_ID_list)

# read the RealDiseaseCodesForPatients.txt file
ukb_phenotypes_df = pd.read_csv('/net/pr2/projects/plgrid/plggneuromol/matzieb/projects/ifpan-michkor-gr/data/UKB-phenotypes/ukb-phenotypes-merged.tsv', sep='\t')

ukb_phenotypes_df.head()


columns_41270 = [col for col in ukb_phenotypes_df.columns if '41270' in col]


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

icd10_codes_41270 = get_phenotype_data_for_field(ukb_phenotypes_df, field_id="41270")

# save icd10_codes_41270 to the file
icd10_codes_41270.to_csv('/net/pr2/projects/plgrid/plggneuromol/matzieb/projects/ifpan-michkor-gr/data/UKB-phenotypes/icd10_codes_41270.tsv', sep='\t', index=False)



# get columns names ukb_phenottypes, split by the '.' and get the first element, then count the values
count_field_id = ukb_phenotypes_df.columns.str.split('.').str[0].value_counts()

# filter elements with count > 1, and get field_id, and drop '41270', 41262 field_id
field_ids = count_field_id[count_field_id > 1].index.tolist()
for item in ['41270', '41262', '41280', '41201', '40002', '41202', '41204', '40006', '40001']:
    if item in field_ids:
        field_ids.remove(item)

len(field_ids)

# prepare empty dataframe with field_id and data type
df = pd.DataFrame(columns=['field_id', 'data_type'])
df = pd.DataFrame()


for field_id in field_ids:
    tmp = get_phenotype_data_for_field(ukb_phenotypes_df, field_id=field_id)
    print("Field ID:", field_id)
    print("Shape:", tmp.shape)
    print("Head:", tmp.head())  
    
    columns_field_id = tmp.columns.tolist()
    
    print("Unique values:", tmp[columns_field_id[0]].unique())
    
    unique_element_list = tmp[columns_field_id[0]].unique().tolist()
    # assess type of data
    print("Data type:", tmp[columns_field_id[0]].dtype)
    
    # add to the dataframe, field_id and type of data, data_type save as string
    df = df._append({'field_id': field_id, 
                     'data_type': str(tmp[columns_field_id[0]].dtype),
                     'number_uniq_values': len(unique_element_list),}, 
                    ignore_index=True)
    
    print("")
    
    
# filter the dataframe with the number of unique values
df[df['number_uniq_values'] == 5].shape

# filter string data type
df[df['data_type'] == 'object'].shape

df[df['data_type'] == 'object'][df['number_uniq_values'] > 10][df['number_uniq_values'] < 50]["field_id"].tolist()


# filter int data type
df[df['data_type'] != 'object'].shape

# select intersection of the field_ids and ukb_phenotype_field_ID_list
df[df['data_type'] != 'object']["field_id"].tolist()

[field_ids[i] for i in ukb_phenotype_field_ID_list]

for field_id in list(set(field_ids) & set(ukb_phenotype_field_ID_list)):
    tmp = get_phenotype_data_for_field(ukb_phenotypes_df, field_id=field_id)
    print("Field ID:", field_id)
    print("Shape:", tmp.shape)
    print("Head:", tmp.head())  
    columns_field_id = tmp.columns.tolist()
    print("Unique values:", tmp[columns_field_id[0]].unique())
    # assess type of data
    print("Data type:", tmp[columns_field_id[0]].dtype)
    print("")
    

    

################################################################################
# convert df to one hot encoding for each column, from each column get unique values and create new columns with the unique values, without Nan values,
# then iterate over the unique values and create new columns with the prefix of the column name and the unique value,
# then use np.where to check if the value in the column is equal to the unique value, if yes, fill the column with 1, if not fill with 0
# then drop the original column with the unique values
def one_hot_encoding_for_df(df):
    for col in df.columns:
        unique_values = df[col].unique()
        for value in unique_values:
            new_col_name = col + "_" + str(value)
            df[new_col_name] = np.where(df[col] == value, 1, 0)
        df.drop(col, axis=1, inplace=True)



# and fill the columns with 1 or 0



def one_hot_encoding_for_df(df):
    for col in df.columns:
        unique_values = df[col].unique()
        print(unique_values)
        for val in unique_values:
            new_col_name = col + "_" + str(val)
            df[new_col_name] = np.where(df[col] == val, 1, 0)
        df.drop(col, axis=1, inplace=True)
    return df


one_hot_encoding_for_df(df = tmp)

ukb_phenotypes_df.columns.tolist()





