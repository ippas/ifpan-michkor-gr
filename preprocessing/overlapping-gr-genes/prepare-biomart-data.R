if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
  
BiocManager::install("biomaRt")

require(readr)

# require(biomaRt)

#### Step 1: Obtain and combine gene data from BioMart ####

# to analysis was used gene from: Ensembl Genes 110

# The BioMart datasets for the three animals
datasets <- list(
  human = "hsapiens_gene_ensembl",
  mouse = "mmusculus_gene_ensembl",
  rat = "rnorvegicus_gene_ensembl"
)

# Define a function to download and process the data from BioMart
download_biomart_data <- function(dataset, animal) {
  # Connect to the Ensembl BioMart database
  mart <- biomaRt::useMart("ensembl")
  # Select the dataset
  genes <- biomaRt::useDataset(dataset, mart)
  # Query the data
  data <- biomaRt::getBM(attributes = c('external_gene_name', 'gene_biotype'), mart = genes)
  # Rename and add animal column
  data <- data %>%
    rename(gene_name = external_gene_name) %>%
    mutate(animal = animal)
  
  return(data)
}

# Apply the function to each dataset and bind the results together
biomart_genes <- purrr::map2_dfr(datasets, names(datasets), download_biomart_data)

# Filter the biomart_genes for protein coding genes
biomart_genes_protein <- biomart_genes %>% 
  filter(gene_biotype == "protein_coding")

# check number of genes for each animal
biomart_genes_protein %>% 
  group_by(animal) %>% 
  nest() %>% 
  mutate(number_genes = map(data, ~nrow(.x))) %>% 
  unnest(number_genes)

#### Step 2: Obtain and combine gene data from Gene Ontology ####

# The URLs for the Gene Ontology data and corresponding animal names
urls <- c(
  "http://geneontology.org/gene-associations/goa_human.gaf.gz",
  "http://current.geneontology.org/annotations/mgi.gaf.gz",
  "http://current.geneontology.org/annotations/rgd.gaf.gz"
)
animals <- c("human", "mouse","rat")

# Define a function to download and process the Gene Ontology data
download_gene_ontology_data <- function(url, animal) {
  # Read the data
  data <- readr::read_tsv(url, skip = 1, comment = "!") 
  # Select the third column (gene names), remove duplicates, and add animal column
  data <- data %>%
    select(3) %>%
    unique() %>%
    rename(gene_name = `...3`) %>%
    mutate(animal = animal)
  
  return(data)
}

# Apply the function to each URL and bind the results together
gene_ontology_genes <- purrr::map2_df(urls, animals, download_gene_ontology_data)

# check number of genes for each animal
gene_ontology_genes %>% 
  group_by(animal) %>% 
  nest() %>% 
  mutate(number_genes = map(data, ~nrow(.x))) %>% 
  unnest(number_genes)



#### Step 3: Filter the BioMart genes based on the Gene Ontology data ####

# Define a function to filter the BioMart data based on the Gene Ontology data
filter_genes <- function(animal) {
  # Select the rows for the current animal and check if gene_name is in the Gene Ontology data
  data <- biomart_genes_protein %>%
    filter(animal == !!animal) %>%
    filter(gene_name %in% gene_ontology_genes[gene_ontology_genes$animal == !!animal, ]$gene_name)
  
  return(data)
}

# Apply the function to each animal and bind the results together
biomart_genes_in_go <- purrr::map_dfr(animals, filter_genes)

# check number of genes for each animal
biomart_genes_in_go %>% 
  group_by(animal) %>% 
  nest() %>% 
  mutate(number_genes = map(data, ~nrow(.x))) %>% 
  unnest(number_genes)

