from functions_dictionary import *

## LOAD START FILES
source_file = r'split_file_40.tsv'
src = 'nervous_290180230_40'
mgi_file_path = r'..\MGI_EntrezGene.tsv'


#SOURCE FILE
publication_df = pd.read_csv(source_file, sep='\t', keep_default_na=False, usecols=[1,2,3,6], skiprows=1, names=['pmid','gene_name','ensembl_id','species'])
publication_geneName = publication_df['gene_name'].tolist()
ensembl_gene = publication_df['ensembl_id'].tolist()
source_number = str(publication_df['pmid'][1])
species = str(publication_df['species'][6])
print(species)


# M`GI
mgi_df = pd.read_csv(mgi_file_path, sep='\t', header=None, usecols=[0, 1, 9], skiprows=1, names=['GeneName', 'EntrezGene', 'Alias'])
ls_row_1 = mgi_df['GeneName'].tolist()
ls_row_2 = mgi_df['EntrezGene'].tolist()
ls_row_10 = mgi_df['Alias'].tolist()


# DATASET
## biomart_server
server = BiomartServer('http://www.ensembl.org/biomart')
if species == 'mouse':
    dataset = server.datasets['mmusculus_gene_ensembl']
    obiekt = musMusculus()
elif species == 'rat':
    dataset = server.datasets['rnorvegicus_gene_ensembl']
    obiekt = rattusNorvegicus()
elif species == 'human':
    dataset = server.datasets['hsapiens_gene_ensembl']
    obiekt = homoSapiens()
    print(dataset)
else:
    print('ERROR: animal not found')
    


# CHANGE_VAR
gene_list_number = 801
source = 'pmid:' + source_number        # => cluster   
gene_list_id = 'all_significant_genes_' + source_number


# LISTS
dictionary = {}
ls_geneDictionaries = []
ls_notResponse = []
ls_notResponse_after = []


alias = 'NA'
info = 'NA'
previous_gene_name = None
previous_gene_data = {}


for i in range(len(publication_geneName)):
    gene_name = publication_geneName[i]
    
    # VARIABLES
    index = i + 1
    ensembl_gene_id_temp = []
    ensembl_transcript_id_temp = []
    refseq_mrna_temp = []
    ensembl_gene_id = ensembl_gene[i]
    print(ensembl_gene_id)
    
    if gene_name == '':
        continue
    
    elif gene_name == 'NA':
        ensembl_gene_id = str('NA')
        ensembl_transcript_id = str('NA')
        refseq_mrna = str('NA')
        hgnc_symbol = str('NA')
        temp_gene_dictionary = gene_dictionary(index,
                                            gene_name,
                                            gene_list_number,
                                            gene_list_id,
                                            source,
                                            ensembl_gene_id,
                                            ensembl_transcript_id, 
                                            refseq_mrna,
                                            hgnc_symbol,
                                            alias,
                                            info)
        ls_geneDictionaries.append(temp_gene_dictionary)
        print(temp_gene_dictionary)
        
    elif gene_name == previous_gene_name:
        # Use the same data as the previous gene
        ensembl_gene_id = previous_gene_data.get('ensembl_gene_id', '')
        ensembl_transcript_id = previous_gene_data.get('ensembl_transcript_id', '')
        refseq_mrna = previous_gene_data.get('refseq_mrna', '')
        hgnc_symbol = previous_gene_data.get('hgnc_symbol', '')
        
        temp_gene_dictionary = gene_dictionary(index,
                                            gene_name,
                                            gene_list_number,
                                            gene_list_id,
                                            source,
                                            ensembl_gene_id,
                                            ensembl_transcript_id, 
                                            refseq_mrna,
                                            hgnc_symbol,
                                            alias,
                                            info)
        ls_geneDictionaries.append(temp_gene_dictionary)
        print(temp_gene_dictionary)
        
        
    else:
        ls_biomartParameters = obiekt.biomartParameters(gene_name, dataset)
        if not ls_biomartParameters:
            ls_notResponse.append(gene_name)
            print(ls_notResponse)


        for ls in ls_biomartParameters:
            for j in range(3):
                if len(ls) < (j+1):
                    ls.append('NA')
            ensembl_gene_id_temp.append(ls[0])
            ensembl_transcript_id_temp.append(ls[1])
            refseq_mrna_temp.append(ls[2])  

        ensembl_gene_id = '|'.join(list(set(filter(None, ensembl_gene_id_temp))))
        ensembl_transcript_id = '|'.join(list(set(filter(None, ensembl_transcript_id_temp))))
        refseq_mrna = '|'.join(list(set(filter(None, refseq_mrna_temp))))
        
        
        # ORTHOLOGS
        if species == 'human':
            hgnc_symbol = gene_name
        
        else:
            ls_biomartHumanOrthologs = obiekt.biomartHumanOrthologs(gene_name, dataset)
            if not ls_biomartHumanOrthologs:
                hgnc_symbol = 'NA'
            else:
                hgnc_symbol = ls_biomartHumanOrthologs[0][0]
            
        
        temp_gene_dictionary = gene_dictionary(index,
                                            gene_name,
                                            gene_list_number,
                                            gene_list_id,
                                            source,
                                            ensembl_gene_id,
                                            ensembl_transcript_id, 
                                            refseq_mrna,
                                            hgnc_symbol,
                                            alias,
                                            info)

        ls_geneDictionaries.append(temp_gene_dictionary)
        print(temp_gene_dictionary)
        
        previous_gene_data = {
            'ensembl_gene_id': ensembl_gene_id,
            'ensembl_transcript_id': ensembl_transcript_id,
            'refseq_mrna': refseq_mrna,
            'hgnc_symbol': hgnc_symbol
        }
    
    previous_gene_name = gene_name

### SCORES 

# from DICTIONARY to file
data = ls_geneDictionaries
dictionary = r'.\Dictionary_' + src + '.tsv'
with open(dictionary, 'w') as file:
    for line in data:
        file.write(str(line) + '\n')

# from ls_notResponse to file
data1 = ls_notResponse
temp_notResponse = r'.\notResponse_' + src + '.tsv'
with open(temp_notResponse, 'w') as file:
    for line in data1:
        file.write(str(line) + '\n')
        print(line)

# LOAD notRESPONSE
ls_notResponse_v2 = pd.read_csv(temp_notResponse, sep='\t', header=None)[0].tolist()
data2 = alias_and_official(ls_notResponse_v2, ls_row_10, ls_row_2, ls_row_1,ls_notResponse_after)
df2 = pd.DataFrame(data2)
print(df2)

# from alias_and_official to file
temp_responseMgi = r'.\responsewithMGI_' + src + '.tsv'
df2.to_csv(temp_responseMgi, sep="\t", index=False, header=None)

# LOAD responsewithMGI
df_response_with_mgi = pd.read_csv(temp_responseMgi, sep='\t')
mgi_id = df_response_with_mgi['mgi_id'].tolist()
gene_name_v2 = df_response_with_mgi['gene_name'].tolist()



### SECOND DICTIONARY
for i in range(len(mgi_id)):
    mgi_ID = mgi_id[i]
    gene_name = gene_name_v2[i]
    
    if gene_name != '':

        # variables
        ensembl_gene_id_temp = []
        ensembl_transcript_id_temp = []
        refseq_mrna_temp = []
        

        ls_biomartParameters_mgi = obiekt.biomartParameters_mgi(mgi_ID, dataset)
        
        if not ls_biomartParameters_mgi:
            print(gene_name)
            ls_notResponse_after.append(gene_name)

        for ls in ls_biomartParameters_mgi:
            for j in range(3):
                if len(ls) < (j+1):
                    ls.append('NA')
            ensembl_gene_id_temp.append(ls[0])
            ensembl_transcript_id_temp.append(ls[1])
            refseq_mrna_temp.append(ls[2])  

        ensembl_gene_id = '|'.join(list(set(filter(None, ensembl_gene_id_temp))))
        ensembl_transcript_id = '|'.join(list(set(filter(None, ensembl_transcript_id_temp))))
        refseq_mrna = '|'.join(list(set(filter(None, refseq_mrna_temp))))
        
        if species == 'human':
            hgnc_symbol = gene_name
        else:
            ls_biomartHumanOrthologs = obiekt.biomartHumanOrthologs_mgi(mgi_ID, dataset)
            if not ls_biomartHumanOrthologs:
                hgnc_symbol = 'NA'
            else:
                hgnc_symbol = ls_biomartHumanOrthologs[0][0]       
        
        
        for i, line in enumerate(ls_geneDictionaries):
            columns = line.split('\t')
            desired_gene_name = columns[1].lower()
            gene_name_lower = str(gene_name).lower()
            
            if desired_gene_name == gene_name_lower:
                print(desired_gene_name, gene_name_lower)
                columns[5] = ensembl_gene_id
                columns[6] = ensembl_transcript_id
                columns[7] = refseq_mrna
                columns[9] = hgnc_symbol
            else:
                ls_notResponse_after.append(gene_name)

# Join the modified columns back line
            ls_geneDictionaries[i] = '\t'.join(columns)

ls_notResponse_after = list(set(ls_notResponse_after))


### SCORES 

# from SECOND DICTIONARY to file
df3 = ls_geneDictionaries
headers_line = ['index','gene_name','gene_list_number','gene_list_index','source','ensembl_gene_id','ensembl_transcript_id','refseq_mrna_id','hgnc_symbol','alias','info']
temp_secDictionary = r'.\secondDictionary_' + src + '.tsv'
with open(temp_secDictionary, 'w') as file:
    file.write('\t'.join(headers_line) + '\n')
    for line in df3:
        file.write(str(line) + '\n')

# from ls_notResponse_after to file
df4 = ls_notResponse_after
temp_lsnotResponseafter = r'.\ls_notResponse_after_' + src + '.tsv'
with open(temp_lsnotResponseafter, 'w') as file:
    for line in df4:
        file.write(str(line) + '\n')

## Add ALIAS
dictionary_file_path = temp_secDictionary
alias_file_path = r'.\withAlias_' + src + '.tsv'
updateCellswithAlias(mgi_file_path, dictionary_file_path, alias_file_path)

# Add INFO
info_file_path = r'.\withINFO_' + src + '.tsv'
updateCellswithINFO(source_file, alias_file_path, info_file_path)
