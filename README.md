# GRTraits : Database of GR regulated genes

#### Project logline (technique, organism, tissue type)

Identification of GR-dependent genes
Organisms: mouse, human, rat

## Methods

This sections should be a description of preprocessin and analysis ready to be included in the publication

## Preprocessing

Details of file preprocessing

FILE DICTIONARY_:
-> sample script adapted to the set of publications

IMPORT
block used to import the necessary packages and modules

DATASET
Biomart server data set for connection
    Human: 'hsapiens_gene_ensembl'
    Mouse: 'mmusculus_gene_ensembl'
    Rat: 'rnorvegicus_gene_ensembl'

FUNCTIONS - gene_dictionary - creates data dictionaries - biomartParameters - query using gene_name value (mgi_symbol) - biomartHumanOrthologs - separate query for HGNC because you cannot get data from two different attribute cards in one query - biomartParametersbyEnsembl - query using ensembl_id value when gene_name is "NA" - biomartHumanOrthologsbyEnsembl - separate query for HGNC by Ensembl - alias_and_official - function which is loooking for alias and mgi_id for our gene_name value in
 MGI_Entrez database if response from Biomart was empty - biomartParameters_synonim - query using gene_synonim value - biomartHumanOrthologs_synonim - separate query for HGNC using gene_synonim value - updateCellswithAlias - function which updates alias column using gtf file for each organism in created dictionaries - updateCellswithINFO - function which updates info column using publication source file.
It compares created dictionaries and publication source file by gene name and add information form different source file columns

LOAD FILES
load workbook, specific worksheet and variables for next steps

DICTIONARY
part of code that completes variables and queries to create dictionaries of genes

SCORES
saving partial scores

SECOND DICTIONARY
completing dictionaries with values from MGI_Entrez database found by functiom alias_and_official

Add ALIAS
use of updateCellswithAlias function

Add INFO
use of updateCellswithINFO function

FILE change_xlsx_TSV:
-> format conversion

## Analysis

Details of analysis

File DICTIONARY_nervous-cells is adapted to search information of rat and add information from specyfic publication source file. To use this script for mouse or human is needed to once configure some information inside because there are many differences between the publications.

Variables: to load specific variables, it's necessary to adapt section:

DICTIONARY SECTION
'gene_list_number', 'gene_list_id', 'source', 'organism'
DATASET SECTION
'dataset' for:
    Human: 'hsapiens_gene_ensembl'
    Mouse: 'mmusculus_gene_ensembl'
    Rat: 'rnorvegicus_gene_ensembl'

Files: to load specific source file, it's necessary to adapt section:

PUBLICATION SOURCE
wb (workbook) and ws (worksheet) variables

Add INFO
source_file_path and source_sheet variables

For Human organism also is needed to turn off functions related to searchig HGNC, and add variable hgnc = gene_name in two places (variables at DICTIONARY and SECOND DICTIONARY section)

_notes: all files included in the repo need to be referenced, either in README or other .md files. The analysis has to be fully reproducible, in principle the repo should contain code + description of how to run it while data and results kept outside_

## About this template

Directories:

- _root_ - README.md, \*.Rproj, general configuration files, etc.
- raw - raw data
- preprocessing - scripts
- data - useful data, created by scripts/tools/preprocessing
- analysis - analysis source code
- results - output ready to present

## Raw

Dataframe:

GTF Files:
https://www.ncbi.nlm.nih.gov/datasets/genome/

Rat:
[mRatBN7.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_015227675.2/)
Mouse:
[GRCm39](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001635.27/)
Human:
[GRCh38.p14](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/)

- https://www.omicspred.org/Scores/Metabolon/INTERVAL
- https://www.omicspred.org/Scores/Nightingale/INTERVAL
- ifpan-GR-database-papers.xlsx (matzieb)
- 12868_2017_352_MOESM1_ESM_02.xlsx (michkor publication)
- Clusters\_-_Symbols (marpiech publication)
