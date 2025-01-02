from functions_dictionary import *

## LOAD START FILES
source_file = r'bronchial-28116096.tsv'
src = 'bronchial-28116096'
gene_list_number = 922
#SOURCE FILE
publication_df = pd.read_csv(source_file, sep='\t', keep_default_na=False, usecols=[1,2,3,6], skiprows=1, names=['pmid','gene_name','ensembl_id','species'])
publication_geneName = publication_df['gene_name'].tolist()
ensembl_gene = publication_df['ensembl_id'].tolist()
source_number = str(publication_df['pmid'][1])
species = str(publication_df['species'][2])
print(species)



# DATASET
## biomart_server
server = BiomartServer('http://www.ensembl.org/biomart')
if species == 'mouse':
    dataset = server.datasets['mmusculus_gene_ensembl']
    mgi_file_path = r'..\alias_genomes\mouse.gtf'
    obiekt = musMusculus()
elif species == 'rat':
    dataset = server.datasets['rnorvegicus_gene_ensembl']
    mgi_file_path = r'..\alias_genomes\rat.gtf'
    obiekt = rattusNorvegicus()
elif species == 'human':
    dataset = server.datasets['hsapiens_gene_ensembl']
    mgi_file_path = r'..\alias_genomes\human.gtf'
    obiekt = homoSapiens()
    print(dataset)
else:
    print('ERROR: animal not found')
    


# CHANGE_VAR
source = 'pmid:' + source_number        # => cluster   
gene_list_id = 'all_significant_genes_' + source_number


# LISTS
dictionary = {}
ls_geneDictionaries = []
ls_notResponse = []
ls_notResponse_after = []


alias = 'NA'
info = 'NA'

for i in range(len(publication_geneName)):
        
    # VARIABLES
    gene_name = publication_geneName[i]
    print(gene_name)
    index = i + 1
    external_gene_name_temp = []
    ensembl_gene_id_temp = []
    ensembl_transcript_id_temp = []
    refseq_mrna_temp = []
    ensembl_gene_id = ensembl_gene[i]
    print(ensembl_gene_id)
    
    if gene_name == '':
        continue
    elif gene_name == 'NA':
        ls_biomartParametersbyEnsembl = obiekt.biomartParametersbyEnsembl(ensembl_gene_id, dataset)
        
        
        for ls in ls_biomartParametersbyEnsembl:
            for j in range(3):
                if len(ls) < (j+1):
                    ls.append('NA')
            external_gene_name_temp.append(ls[0])
            ensembl_transcript_id_temp.append(ls[1])
            refseq_mrna_temp.append(ls[2])  

        external_gene_name = '|'.join(list(set(filter(None, external_gene_name_temp))))
        ensembl_transcript_id = '|'.join(list(set(filter(None, ensembl_transcript_id_temp))))
        refseq_mrna = '|'.join(list(set(filter(None, refseq_mrna_temp))))
        
        hgnc_symbol = external_gene_name
        temp_gene_dictionary = gene_dictionary(index,
                                            external_gene_name,
                                            gene_list_id,
                                            gene_list_number,
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
                                            gene_list_id,
                                            gene_list_number,
                                            source,
                                            ensembl_gene_id,
                                            ensembl_transcript_id, 
                                            refseq_mrna,
                                            hgnc_symbol,
                                            alias,
                                            info)

        ls_geneDictionaries.append(temp_gene_dictionary)
        print(temp_gene_dictionary)

### SCORES 

# from DICTIONARY to file
data = ls_geneDictionaries
dictionary = r'.\Dictionary\Dictionary_' + src + '.tsv'
with open(dictionary, 'w') as file:
    for line in data:
        file.write(str(line) + '\n')

# from ls_notResponse to file
data1 = ls_notResponse
temp_notResponse = r'.\ls_notResponse\notResponse_' + src + '.tsv'

with open(temp_notResponse, 'w') as file:
    for symbol in data1:
        file.write(str(symbol) + '\n')
        print(symbol)

if ls_notResponse:
# LOAD notRESPONSE
    ls_notResponse_v2 = pd.read_csv(temp_notResponse, sep='\t', header=None)[0].tolist()
    ls_row_1, ls_row_10 = obiekt.get_gene_lists()
    data2 = alias_and_official(ls_notResponse_v2, ls_row_10, ls_row_1)
    df2 = pd.DataFrame(data2)
    print(df2)
else:
    df2 = pd.DataFrame()

# from alias_and_official to file
temp_responseSynonim = r'.\response_withSynonim\responsewithSynonim_' + src + '.tsv'
df2.to_csv(temp_responseSynonim, sep="\t", index=False, header=None)

# LOAD responsewithSynonim
if ls_notResponse:
    df_temp_responseSynonim = pd.read_csv(temp_responseSynonim, sep='\t')
    synonim_id = df_temp_responseSynonim['gene_synonim'].tolist()
    gene_name_v2 = df_temp_responseSynonim['gene_name'].tolist()
else:
    synonim_id = []
    gene_name_v2 = []


### SECOND DICTIONARY
for i in range(len(synonim_id)):
    synonim_ID = synonim_id[i]
    gene_name = gene_name_v2[i]
    
    if gene_name != '':

        # variables
        ensembl_gene_id_temp = []
        ensembl_transcript_id_temp = []
        refseq_mrna_temp = []
        

        ls_biomartParameters_synonim = obiekt.biomartParameters_synonim(synonim_ID, dataset)
        
        if not ls_biomartParameters_synonim:
            print(gene_name)
            ls_notResponse_after.append(gene_name)

        for ls in ls_biomartParameters_synonim:
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
            ls_biomartHumanOrthologs = obiekt.biomartHumanOrthologs_synonim(synonim_ID, dataset)
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


            ls_geneDictionaries[i] = '\t'.join(columns)

ls_notResponse_after = list(set(ls_notResponse_after))


### SCORES 

# from SECOND DICTIONARY to file
df3 = ls_geneDictionaries
headers_line = ['index','gene_name','gene_list_index','gene_list_number','source','ensembl_gene_id','ensembl_transcript_id','refseq_mrna_id','hgnc_symbol','alias','info']
temp_secDictionary = r'.\secondDictionary\secondDictionary_' + src + '.tsv'
with open(temp_secDictionary, 'w') as file:
    file.write('\t'.join(headers_line) + '\n')
    for line in df3:
        file.write(str(line) + '\n')

# from ls_notResponse_after to file
df4 = ls_notResponse_after
temp_lsnotResponseafter = r'.\notResponse_after\ls_notResponse_after_' + src + '.tsv'
with open(temp_lsnotResponseafter, 'w') as file:
    for line in df4:
        file.write(str(line) + '\n')

## Add ALIAS
dictionary_file_path = temp_secDictionary
alias_file_path = r'.\withAlias\withAlias_' + src + '.tsv'
updateCellswithAlias(mgi_file_path, dictionary_file_path, alias_file_path)

# Add INFO
info_file_path = r'.\withINFO\withINFO_' + src + '.tsv'
updateCellswithINFO(source_file, alias_file_path, info_file_path)

# Fix empty

with open(info_file_path, 'r', newline='') as file:
    rows = [[col if col else 'NA' for col in row] for row in csv.reader(file, delimiter='\t')]

with open(info_file_path, 'w', newline='') as file:
    csv.writer(file, delimiter='\t').writerows(rows)


print('finish')