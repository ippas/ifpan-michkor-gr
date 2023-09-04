# IMPORTS
import csv
from biomart import BiomartServer
import pandas as pd
import numpy as np



# FUNCTIONS
def gene_dictionary(index,
                 gene_name,
                 gene_list_number,
                 gene_list_id,
                 source,
                 ensembl_gene_id,
                 ensembl_transcript_id,
                 refseq_mrna,
                 hgnc_symbol,
                 alias,
                 info):

    dictionary = '\t'.join([
        str(index) if index else 'NA',
        str(gene_name) if gene_name else 'NA',
        str(gene_list_number) if gene_list_number else 'NA',
        str(gene_list_id) if gene_list_id else 'NA',
        str(source) if source else 'NA',
        str(ensembl_gene_id) if ensembl_gene_id else 'NA',
        str(ensembl_transcript_id) if ensembl_transcript_id else 'NA',
        str(refseq_mrna) if refseq_mrna else 'NA',
        str(hgnc_symbol) if hgnc_symbol else 'NA',
        str(alias) if alias else 'NA',
        str(info) if info else 'NA'
    ])


    return dictionary

  
def biomartParameters(mgi_symbol, dataset):
    attributes = ['ensembl_gene_id',
                  'ensembl_transcript_id',
                  'refseq_mrna']
    filters = {'mgi_symbol':[mgi_symbol]}            # gene_name = mgi_symbol
    response = dataset.search({'attributes':attributes,'filters':filters})
    
    # response convertion
    values = [line.split("\t") for line in response.text.split("\n") if line.strip()]
        
    return values 


def biomartHumanOrthologs(mgi_symbol, dataset):
    attributes = ['hsapiens_homolog_associated_gene_name']
    filters = {'mgi_symbol':[mgi_symbol]}
    response = dataset.search({'attributes':attributes,'filters':filters})
    
    # response_convertion
    values = [line.split("\t") for line in response.text.split("\n") if line.strip()]
        
    return values 


def alias_and_official(ls_notResponse,ls_row_10,ls_row_2,ls_row_1, ls_notResponse_after):

    ls_response = []
    ls_response_2 = []
    ls_response_3 = []
    
    for notResponse in ls_notResponse:
        word = str(notResponse).lower()
        for i in range(len(ls_row_10)):
            if ls_row_10[i]:
                temp = str(ls_row_10[i]).lower()
                t_strings = temp.split("|")
                if word in t_strings:
                    ls_response.append([word, i, ls_row_1[i]])

    for notResponse_2 in ls_notResponse:
        word = str(notResponse_2).lower()
        for i in range(len(ls_row_2)):
            if ls_row_2[i]:
                temp = str(ls_row_2[i]).lower()
                if word == temp:
                    ls_response_2.append([word, i, ls_row_1[i]])
                    
    
    ls_response.insert(0, ['gene_name','official/alias_index','mgi_id'])

    ls_response_3 = [ls_response[0]]  # Initialize with the header from ls_response

    gene_names_seen = set()

    for el in ls_response_2[0:]:
        gene_name = el[0]
        if gene_name not in gene_names_seen:
            gene_names_seen.add(gene_name)
            ls_response_3.append(el)

    for el in ls_response[1:]:
        gene_name = el[0]
        if gene_name not in gene_names_seen:
            gene_names_seen.add(gene_name)
            ls_response_3.append(el)
      
    return ls_response_3

server = BiomartServer('http://www.ensembl.org/biomart')       
dataset = server.datasets['mmusculus_gene_ensembl']

def biomartParameters_mgi(mgi_id, dataset):

    attributes = ['ensembl_gene_id',
                  'ensembl_transcript_id',
                  'refseq_mrna']
    filters = {'mgi_id':[mgi_id]}                     # gene_name = mgi_id
    response = dataset.search({'attributes':attributes,'filters':filters})

    # response_convertion
    values = [line.split("\t") for line in response.text.split("\n") if line.strip()]
        
    return values 


def biomartHumanOrthologs_mgi(mgi_id, dataset):
    attributes = ['hsapiens_homolog_associated_gene_name']
    
    filters = {'mgi_id':[mgi_id]}
    response = dataset.search({'attributes':attributes,'filters':filters})
    
    # response_convertion
    values = [line.split("\t") for line in response.text.split("\n") if line.strip()]
        
    return values 

'''
def get_regulation_for_gene(gene_name, source_file):
    with open(source_file, 'r') as f:
        source_data = [line.strip().split('\t') for line in f.readlines()]

    for row in source_data[1:]:
        if row[3] == gene_name:                 # Compare FirstOfName to gene_name
            fold_change = float(row[log2Ratio])        # fold change is in the 34th column
            if fold_change < 0:
                return 'DOWN' 
            else:
                return 'UP'
'''


def updateCellswithAlias(mgi_file_path, dictionary_file_path, alias_file_path):
    wb_mgi = pd.read_csv(mgi_file_path, sep='\t', header=None, keep_default_na=False, usecols=[1, 9], skiprows=2, names=['g_name','alias'])
    mgi_dict = dict(zip(wb_mgi['g_name'], wb_mgi['alias']))
    
    wb_dictionary = pd.read_csv(dictionary_file_path, keep_default_na=False, sep='\t', header=0)  # Assuming the first row is the header

    def debug_lambda(x):
        result = mgi_dict.get(x, 'NA')
        if result == '' or None:
            print(f"Debug: x={x}, mgi_dict.get(x)=NaN, setting result to 'NA'")
            result = 'NA'
        else:
            print(f"Debug: x={x}, mgi_dict.get(x)={result}, result={result}")
        return result

    wb_dictionary['alias'] = wb_dictionary['gene_name'].apply(debug_lambda)
    wb_dictionary.to_csv(alias_file_path, sep='\t', index=False)
   

 
def updateCellswithINFO(source_file_path, alias_file_path, info_file_path):
    
    # Read source_file
    with open(source_file_path, 'r') as f:
        csv_reader = csv.reader(f, delimiter='\t')
        headers = next(csv_reader)
        info_data = [row for row in csv_reader]

    ls_row_g_name = [row[headers.index('gene_name')] for row in info_data]

    # Read alias_file
    with open(alias_file_path, 'r') as f:
        csv_reader = csv.reader(f, delimiter='\t')
        alias_headers = next(csv_reader)
        alias_data = [row for row in csv_reader]

    ls_row_2 = [row[1] for row in alias_data] 

    # Comparing both
    for id, alias in enumerate(ls_row_2):
        for j, g_name in enumerate(ls_row_g_name):
            try:
                if g_name.lower() == alias.lower():
                    print(f"Found a match for {g_name.lower()} and {alias.lower()}")
                    new_data = '|'.join([f"{headers[i]}:" + str(info_data[j][i]) for i in range(4, len(headers))])
                    alias_data[id][10] = new_data 

            except Exception as e:
                print(f"An exception occurred: {e}")
                
    # Write to file info
    with open(info_file_path, 'w', newline='') as f:
        csv_writer = csv.writer(f, delimiter='\t')
        csv_writer.writerow(alias_headers)
        for row in alias_data:
            csv_writer.writerow(row)
