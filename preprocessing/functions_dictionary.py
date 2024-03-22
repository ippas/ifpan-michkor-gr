# IMPORTS
import csv
from biomart import BiomartServer
import pandas as pd
import numpy as np
from pybiomart import Dataset


# GENERAL FUNCTIONS 
def gene_dictionary(index,
                    gene_name,
                    gene_list_id,
                    gene_list_number,
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
            str(gene_list_id) if gene_list_id else 'NA',
            str(gene_list_number) if gene_list_number else 'NA',
            str(source) if source else 'NA',
            str(ensembl_gene_id) if ensembl_gene_id else 'NA',
            str(ensembl_transcript_id) if ensembl_transcript_id else 'NA',
            str(refseq_mrna) if refseq_mrna else 'NA',
            str(hgnc_symbol) if hgnc_symbol else 'NA',
            str(alias) if alias else 'NA',
            str(info) if info else 'NA'
        ])

        return dictionary


def alias_and_official(ls_notResponse, ls_row_10, ls_row_1):

    ls_row_1_processed = [(str(item).casefold(), item, i) for i, item in enumerate(ls_row_1) if item]
    ls_row_10_processed = [(str(item).casefold(), item, i) for i, item in enumerate(ls_row_10) if item]

    ls_response_combined = [['gene_name', 'official/alias_index', 'gene_synonim']]
    gene_names_seen = set()

    for notResponse in ls_notResponse:
        word = str(notResponse).casefold()
        matched = False

        for processed_list in (ls_row_1_processed, ls_row_10_processed):
            for temp, original, i in processed_list:
                if word in temp.split(", ") or word == temp:

                    entry = [notResponse, i, original]

                    if entry[0] not in gene_names_seen:
                        gene_names_seen.add(entry[0])
                        ls_response_combined.append(entry)
                        matched = True
                        break 
            if matched:
                break 

    return ls_response_combined


def read_gtf_and_extract_genes_synonyms(mgi_file_path):
    # Read and process the GTF file
    mgi_df = pd.read_csv(mgi_file_path, sep='\t', header=None, comment='#', 
                         names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes'])
    mgi_df['GeneName'] = mgi_df['attributes'].str.extract('gene_id "([^"]+)"')
    mgi_df['GeneSynonym'] = mgi_df['attributes'].str.findall('gene_synonym "([^"]+)"').apply(lambda x: ', '.join(x) if x else '')

    
    gene_synonyms_aggregated = {}
    for _, row in mgi_df.iterrows():
        gene_name = row['GeneName'].strip().casefold()
        synonyms = set(row['GeneSynonym'].split('|')) if row['GeneSynonym'] else set()
        
        if gene_name in gene_synonyms_aggregated:
            gene_synonyms_aggregated[gene_name].update(synonyms)
        else:
            gene_synonyms_aggregated[gene_name] = synonyms

    
    for gene_name, synonyms in gene_synonyms_aggregated.items():
        gene_synonyms_aggregated[gene_name] = '|'.join(synonyms)

    return gene_synonyms_aggregated

def updateCellswithAlias(mgi_file_path, dictionary_file_path, alias_file_path):
    gene_synonyms_aggregated = read_gtf_and_extract_genes_synonyms(mgi_file_path)

    wb_dictionary = pd.read_csv(dictionary_file_path, keep_default_na=False, sep='\t', header=0)

    def update_alias(gene_name):
        gene_name_cf = gene_name.casefold()
        for gene, synonyms in gene_synonyms_aggregated.items():
            if gene_name_cf == gene or gene_name_cf in synonyms.split('|'):
                return synonyms or 'NA'
        return 'NA'

    wb_dictionary['alias'] = wb_dictionary['gene_name'].apply(update_alias)

    wb_dictionary.to_csv(alias_file_path, sep='\t', index=False)




'''
def updateCellswithAlias(mgi_file_path, dictionary_file_path, alias_file_path):
    
    wb_mgi = pd.read_csv(mgi_file_path, sep='\t', header=None, keep_default_na=False, usecols=[1, 9], skiprows=2, names=['g_name','alias'])
    mgi_dict = dict(zip(wb_mgi['gene_name'], wb_mgi['geme']))
    wb_dictionary = pd.read_csv(dictionary_file_path, keep_default_na=False, sep='\t', header=0)

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
'''    
    
def updateCellswithINFO(source_file_path, alias_file_path, info_file_path):
    
    with open(source_file_path, 'r') as f:
        csv_reader = csv.reader(f,delimiter='\t')
        headers = next(csv_reader)
        info_data = [row for row in csv_reader]
        
    with open(alias_file_path, 'r') as f:
        csv_reader = csv.reader(f, delimiter='\t')
        alias_headers = next(csv_reader)
        alias_data = [row for row in csv_reader]

    
    for id, alias_row in enumerate(alias_data):
        try:
            info_row = info_data[id]                
            new_data = '|'.join([f"{headers[i]}:" + str(info_row[i]) for i in range(4, len(headers))])
            
            if len(alias_row) > 10:  
                alias_row[10] = new_data
            else:
                print(f"Row {id} in alias_data has fewer than 11 columns.")
        
        except Exception as e:
            print(f"An exception occurred: {e}")

    with open(info_file_path, 'w', newline='') as f:
        csv_writer = csv.writer(f, delimiter='\t')
        csv_writer.writerow(alias_headers)
        for row in alias_data:
            csv_writer.writerow(row)
            
    
## FUNCTIONS SETS - ONE CLASS FOR SPECIES
class musMusculus:
    
    def __init__(self):
        self.mgi_file_path = r'..\alias_genomes\mouse.gtf'
        self.mgi_df = None

    def read_gtf_file(self):
        self.mgi_df = pd.read_csv(self.mgi_file_path, sep='\t', header=None, comment='#', 
                                names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes'])
        self.mgi_df['GeneName'] = self.mgi_df['attributes'].str.extract('gene_id "([^"]+)"')
        self.mgi_df['GeneSynonym'] = self.mgi_df['attributes'].str.findall('gene_synonym "([^"]+)"').apply(lambda x: ', '.join(x) if x else '')

    def get_gene_lists(self):
        if self.mgi_df is None:
            self.read_gtf_file()

        ls_row_1 = self.mgi_df['GeneName'].tolist()
        ls_row_10 = self.mgi_df['GeneSynonym'].tolist()
        
        return ls_row_1, ls_row_10
    
    # METODY
    def biomartParameters(self, mgi_symbol, dataset):
        attributes = ['ensembl_gene_id',
                    'ensembl_transcript_id',
                    'refseq_mrna']
        filters = {'mgi_symbol':[mgi_symbol]}            # gene_name = mgi_symbol
        response = dataset.search({'attributes':attributes,'filters':filters})
        
        values = [line.split("\t") for line in response.text.split("\n") if line.strip()]
        return values 


    def biomartHumanOrthologs(self, mgi_symbol, dataset):
        attributes = ['hsapiens_homolog_associated_gene_name']
        filters = {'mgi_symbol':[mgi_symbol]}
        response = dataset.search({'attributes':attributes,'filters':filters})
        
        values = [line.split("\t") for line in response.text.split("\n") if line.strip()]
        return values 


    def biomartParameters_synonim(self, external_synonym, dataset):

        attributes = ['ensembl_gene_id',
                    'ensembl_transcript_id',
                    'refseq_mrna']
        filters = {'external_synonym':[external_synonym]}                     # gene_synonim
        response = dataset.search({'attributes':attributes,'filters':filters})

        # response_convertion
        values = [line.split("\t") for line in response.text.split("\n") if line.strip()]
        return values 

    def biomartHumanOrthologs_synonim(self, external_synonym, dataset):
        attributes = ['hsapiens_homolog_associated_gene_name']
        
        filters = {'external_synonym':[external_synonym]}
        response = dataset.search({'attributes':attributes,'filters':filters})
        
        # response_convertion
        values = [line.split("\t") for line in response.text.split("\n") if line.strip()]
        return values
    
    def biomartParametersbyEnsembl(self, ensembl_gene_id, dataset):

        attributes = ['external_gene_name',
                    'ensembl_transcript_id',
                    'refseq_mrna']
        filters = {'ensembl_gene_id':[ensembl_gene_id]}               # gene_name = mgi_symbol
        response = dataset.search({'attributes':attributes,'filters':filters})
        
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

                
class rattusNorvegicus:
    
    def __init__(self):
        self.mgi_file_path = r'..\alias_genomes\rat.gtf'
        self.mgi_df = None
        

    def read_gtf_file(self):

        self.mgi_df = pd.read_csv(self.mgi_file_path, sep='\t', header=None, comment='#', 
                                names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes'])
        self.mgi_df['GeneName'] = self.mgi_df['attributes'].str.extract('gene_id "([^"]+)"')
        self.mgi_df['GeneSynonym'] = self.mgi_df['attributes'].str.findall('gene_synonym "([^"]+)"').apply(lambda x: ', '.join(x) if x else '')

    def get_gene_lists(self):
        if self.mgi_df is None:
            self.read_gtf_file() 

        ls_row_1 = self.mgi_df['GeneName'].tolist()
        ls_row_10 = self.mgi_df['GeneSynonym'].tolist()
        
        return ls_row_1, ls_row_10
    
    def biomartParameters(self, mgi_symbol, dataset):
        attributes = ['ensembl_gene_id',
                    'ensembl_transcript_id',
                    'refseq_mrna']
        filters = {'external_gene_name':[mgi_symbol]}            # gene_name = mgi_symbol
        response = dataset.search({'attributes':attributes,'filters':filters})
        
        values = [line.split("\t") for line in response.text.split("\n") if line.strip()]
        return values 


    def biomartHumanOrthologs(self, mgi_symbol, dataset):
        attributes = ['hsapiens_homolog_associated_gene_name']
        filters = {'external_gene_name':[mgi_symbol]}
        response = dataset.search({'attributes':attributes,'filters':filters})
        
        values = [line.split("\t") for line in response.text.split("\n") if line.strip()]
        return values 

    
    def biomartParametersbyEnsembl(self, ensembl_gene_id, dataset):

        attributes = ['external_gene_name',
                    'ensembl_transcript_id',
                    'refseq_mrna']
        filters = {'ensembl_gene_id':[ensembl_gene_id]}               # gene_name = mgi_symbol
        response = dataset.search({'attributes':attributes,'filters':filters})
        
        values = [line.split("\t") for line in response.text.split("\n") if line.strip()] 
        return values


    def biomartHumanOrthologsbyEnsembl(self, dataset, ensembl_gene_id):

        attributes = ['hsapiens_homolog_associated_gene_name']
        filters = {'ensembl_gene_id':[ensembl_gene_id]}               # gene_name = mgi_symbol
        response = dataset.search({'attributes':attributes,'filters':filters})
        
        values = [line.split("\t") for line in response.text.split("\n") if line.strip()] 
        return values 

    def biomartParameters_synonim(self, external_synonym, dataset):

        attributes = ['ensembl_gene_id',
                    'ensembl_transcript_id',
                    'refseq_mrna']
        filters = {'external_synonym':[external_synonym]}                     # gene_synonim
        response = dataset.search({'attributes':attributes,'filters':filters})

        values = [line.split("\t") for line in response.text.split("\n") if line.strip()]
        return values 

    def biomartHumanOrthologs_synonim(self, external_synonym, dataset):
        attributes = ['hsapiens_homolog_associated_gene_name']
        
        filters = {'external_synonym':[external_synonym]}
        response = dataset.search({'attributes':attributes,'filters':filters})
        
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
    
    
class homoSapiens:
    
    def __init__(self):
        # Initialize the file path upon creating an instance of the class
        self.mgi_file_path = r'..\alias_genomes\human.gtf'
        self.mgi_df = None

    def read_gtf_file(self):
        # Read the file and process it
        self.mgi_df = pd.read_csv(self.mgi_file_path, sep='\t', header=None, comment='#', 
                                  names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes'])
        self.mgi_df['GeneName'] = self.mgi_df['attributes'].str.extract('gene_id "([^"]+)"')
        self.mgi_df['GeneSynonym'] = self.mgi_df['attributes'].str.findall('gene_synonym "([^"]+)"').apply(lambda x: ', '.join(x) if x else '')

    def get_gene_lists(self):
        if self.mgi_df is None:
            self.read_gtf_file()  # Ensure the file is read and processed before attempting to access the data
        
        # Convert to list as required
        ls_row_1 = self.mgi_df['GeneName'].tolist()
        ls_row_10 = self.mgi_df['GeneSynonym'].tolist()
        
        return ls_row_1, ls_row_10
    
    def biomartParameters(self, mgi_symbol, dataset):

        attributes = ['ensembl_gene_id',
                    'ensembl_transcript_id',
                    'refseq_mrna']

        filters = {'external_gene_name':[mgi_symbol]}               # gene_name = mgi_symbol
        response = dataset.search({'attributes':attributes,'filters':filters})
        
        # response_convertion
        values = [line.split("\t") for line in response.text.split("\n") if line.strip()]
        return values   

    def biomartParameters_synonim(self, gene_name, dataset):

        attributes = ['ensembl_gene_id',
                    'ensembl_transcript_id',
                    'refseq_mrna']
        filters = {'uniprot_gn_symbol':[gene_name]}                     # gene_name = mgi_id
        response = dataset.search({'attributes':attributes,'filters':filters})
        
        # response_convertion
        values = [line.split("\t") for line in response.text.split("\n") if line.strip()]
        return values
    
    def biomartParametersbyEnsembl(self, ensembl_gene_id, dataset):

        attributes = ['external_gene_name',
                    'ensembl_transcript_id',
                    'refseq_mrna']
        filters = {'ensembl_gene_id':[ensembl_gene_id]}               # gene_name = mgi_symbol
        response = dataset.search({'attributes':attributes,'filters':filters})
        
        values = [line.split("\t") for line in response.text.split("\n") if line.strip()] 
        return values  
'''

    def biomartHumanOrthologsbyEnsembl(self, ensembl_gene_id, dataset):

        attributes = ['hsapiens_homolog_associated_gene_name']
        filters = {'ensembl_gene_id':[ensembl_gene_id]}               # gene_name = mgi_symbol
        response = dataset.search({'attributes':attributes,'filters':filters})
        
        values = [line.split("\t") for line in response.text.split("\n") if line.strip()] 
        return values 
'''
    
