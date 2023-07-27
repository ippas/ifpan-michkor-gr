if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")

install.packages("circlize")
install.packages("pheatmap")

# load packages
require(tidyverse)
require(magrittr)
require(ComplexHeatmap)
require(circlize)
library(pheatmap)

source("preprocessing/overlapping-gr-genes/gr-gene-database-functions.  ")
getwd()
# read data
gr_gene_database_raw <- read.csv("data/overlapping-gr-genes/geneBase_060723.tsv", sep = "\t") 

gr_gene_database_raw %>% head

gr_gene_database_raw %>%  colnames()

gr_gene_database_raw %>% 
  extract_keys_values("Info", c("tissue", "cell", "enviroment", " treatment", "dose", "time", "Log2Ratio")) -> gr_gene_database

gr_gene_database %>%   
  # tail(100) %>% 
  print_simple_database() %>% 
  mutate(Source = str_replace_all(Source, "PMID: ", "")) %>% 
  filter(Source != "marpiech_tissues_dex") %>% 
  filter(Source != "michkor-cells-dex") %>% 
  mutate(Gene_name = tolower(Gene_name)) %>% 
  mutate(label = paste0(Source, "_", tissue, "_", cell)) -> gr_gene_database_preprocessing


gr_gene_database_preprocessing %>% 
  filter(Gene_name %in% {biomart_genes_in_go %>% .$gene_name %>% tolower() %>% unique()}) %>% 
  filter(label != "34362910_NA_NA")-> gr_gene_filter_biomart
  

split(gr_gene_filter_biomart$Gene_name, gr_gene_filter_biomart$label) -> gr_gene_list
total_genes = biomart_genes_in_go %>% .$gene_name %>% tolower() %>% unique()

# Use the function
results <- perform_chi2_tests(gr_gene_list, total_genes)


results$chi2_value_matrix -> chi2_value_matrix
chi2_value_matrix <- ifelse(chi2_value_matrix > 5000, NA, chi2_value_matrix)
chi2_value_matrix <- ifelse(chi2_value_matrix > 700, 700, chi2_value_matrix)


pheatmap(chi2_value_matrix, cluster_rows = T, cluster_cols = T, display_numbers = T)

