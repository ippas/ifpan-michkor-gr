# List of CRAN packages
cran_packages <- c("tidyverse", "magrittr", "circlize", "pheatmap")

# List of Bioconductor packages
bioc_packages <- c("ComplexHeatmap")

# Function to install and load packages
load_package <- function(pkg, is_bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (is_bioc) {
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
  }
  library(pkg, character.only = TRUE)
}

# Check if BiocManager is installed, if not install
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install and load CRAN packages
lapply(cran_packages, load_package)

# Install and load Bioconductor packages
lapply(bioc_packages, load_package, is_bioc = TRUE)

# Check if biomaRt is installed
if (!requireNamespace("biomaRt", quietly = TRUE)) {
  install.packages("biomaRt")
}

detach(package:biomaRt)

source("preprocessing/overlapping-gr-genes/gr-database-functions.R")
source("preprocessing/overlapping-gr-genes/prepare-biomart-data.R")


# Read data
gr_gene_database_raw <- read.csv("data/overlapping-gr-genes/geneBase_060723.tsv", sep = "\t") 

# Extract specific keys and values from the "Info" column
gr_gene_database_raw %>% 
  extract_keys_values("Info", c("tissue", "cell", "enviroment", "treatment", "dose", "time", "Log2Ratio")) -> gr_gene_database

# Preprocess the gene database
gr_gene_database %>%   
  print_simple_database() %>%
  mutate(Source = str_replace_all(Source, "PMID: ", "")) %>%
  filter(Source != "michkor-cells-dex") %>%
  mutate(Gene_name = tolower(Gene_name)) %>%
  mutate(label = paste0(Source, "_", tissue, "_", cell)) %>%
  mutate(label = ifelse(Source == "marpiech_tissues_dex", paste0(Source, "_", Gene_list_number), label)) %>% 
  filter(Gene_name %in% {filtered_biomart_by_go %>% .$gene_name %>% tolower() %>% unique()}) %>%
  filter(label != "34362910_NA_NA") -> gr_gene_filter_biomart


# Split gene names by label
split(gr_gene_filter_biomart$Gene_name, gr_gene_filter_biomart$label) -> gr_gene_list

# Get unique total genes
total_genes = filtered_biomart_by_go %>% .$gene_name %>% tolower() %>% unique()

# Perform chi-square tests
results <- perform_chi2_tests(gr_gene_list, total_genes)

# Create heatmaps
pheatmap(log2(results$chi2_value_matrix+1), cluster_rows = T, cluster_cols = T, display_numbers = T)
pheatmap(results$p_value_matrix, cluster_rows = T, cluster_cols = T, display_numbers = T)

Heatmap(log2(results$chi2_value_matrix+1), 
        name = "overlap", 
        col = color_palette, 
        show_row_names = TRUE, 
        show_column_names = TRUE)

