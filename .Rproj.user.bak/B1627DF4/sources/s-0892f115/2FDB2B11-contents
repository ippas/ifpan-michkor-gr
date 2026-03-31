#### Step 1: Obtain and Combine Gene Data from BioMart ####
# Utilizing Ensembl Genes 110 for analysis

# Define BioMart datasets for three animals
datasets <- list(
  human = "hsapiens_gene_ensembl",
  mouse = "mmusculus_gene_ensembl",
  rat = "rnorvegicus_gene_ensembl"
)

# Function to download and process data from BioMart
download_biomart_data <- function(dataset, animal) {
  mart <- biomaRt::useMart("ensembl")
  genes <- biomaRt::useDataset(dataset, mart)
  data <- biomaRt::getBM(attributes = c('external_gene_name', 'gene_biotype'), mart = genes)
  data <- data %>%
    rename(gene_name = external_gene_name) %>%
    mutate(animal = animal)
  return(data)
}

# Obtain BioMart genes and filter for protein coding genes
biomart_genes <- purrr::map2_dfr(datasets, names(datasets), download_biomart_data)
biomart_genes_protein <- biomart_genes %>% filter(gene_biotype == "protein_coding")

# Check number of protein coding genes for each animal
biomart_genes_protein %>%
  group_by(animal) %>%
  nest() %>%
  mutate(number_genes = map(data, ~nrow(.x))) %>%
  unnest(number_genes)


#### Step 2: Obtain and Combine Gene Data from Gene Ontology ####
# Define URLs for Gene Ontology data
urls <- c(
  "http://current.geneontology.org/annotations/goa_human.gaf.gz",
  "http://current.geneontology.org/annotations/mgi.gaf.gz",
  "http://current.geneontology.org/annotations/rgd.gaf.gz"
)
animals <- c("human", "mouse", "rat")

# Function to download and process Gene Ontology data
download_gene_ontology_data <- function(url, animal) {
  data <- readr::read_tsv(url, skip = 1, comment = "!")
  data <- data %>%
    select(3) %>%
    unique() %>%
    set_colnames(c("gene_name")) %>%
    mutate(animal = animal)
  return(data)
}

# Combine Gene Ontology data for each animal
gene_ontology_genes <- purrr::map2_df(urls, animals, download_gene_ontology_data)

# Check number of genes for each animal in Gene Ontology data
gene_ontology_genes %>%
  group_by(animal) %>%
  nest() %>%
  mutate(number_genes = map(data, ~nrow(.x))) %>%
  unnest(number_genes)


#### Step 3: Filter BioMart Genes Based on Gene Ontology Data ####
# Function to filter BioMart data based on Gene Ontology data
filter_genes <- function(animal) {
  data <- biomart_genes_protein %>%
    filter(animal == !!animal) %>%
    filter(gene_name %in% gene_ontology_genes[gene_ontology_genes$animal == !!animal, ]$gene_name)
  return(data)
}

# Filter and combine results for each animal
filtered_biomart_by_go <- purrr::map_dfr(animals, filter_genes)

# Check number of filtered genes for each animal
filtered_biomart_by_go %>%
  group_by(animal) %>%
  nest() %>%
  mutate(number_genes = map(data, ~nrow(.x))) %>%
  unnest(number_genes)

# Cleanup: Remove unnecessary variables
rm(datasets, biomart_genes, biomart_genes_protein, urls, animals, gene_ontology_genes)