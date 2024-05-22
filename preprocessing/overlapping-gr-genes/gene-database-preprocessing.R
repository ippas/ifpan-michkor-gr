source("preprocessing/functions/R/install-load-packages.R")

# List of CRAN packages
cran_packages <- c("tidyverse", "magrittr", "circlize", "pheatmap", "yaml", "reshape2", "ComplexHeatmap", "criclize", "RColorBrewer", "grid")

# List of Bioconductor packages
bioc_packages <- c("ComplexHeatmap")

# Function to install and load packages
load_package <- function(pkg, is_bioc = FALSE) {
  # Check if package is already loaded
  if (!requireNamespace(pkg, quietly = TRUE)) {
    # Install from Bioconductor if specified
    if (is_bioc) {
      BiocManager::install(pkg)
    } else { # Install from CRAN
      install.packages(pkg)
    }
  }
  # Load the package
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

# Check if biomaRt is installed, if not install it
if (!requireNamespace("biomaRt", quietly = TRUE)) {
  BiocManager::install("biomaRt")
}

# Detach the biomaRt package
detach(package:biomaRt)

################################################################################
# Read data
# gr_gene_database_raw <- read.csv("data/overlapping-gr-genes/geneBase_060723.tsv", sep = "\t") 
# gr_gene_database_raw <- read.csv("data/overlapping-gr-genes/publikacje_gr_100923.tsv", sep = "\t") 
# gr_gene_database_raw <- read.csv("data/overlapping-gr-genes/publikacje_gr_v5_181023.tsv", sep = "\t")
# gr_gene_database_raw <- read.csv("data/overlapping-gr-genes/gr_geneBase_v8_251023.tsv", sep = "\t") 
# gr_gene_database_raw <- read.csv("data/overlapping-gr-genes/gr_geneBase_v9_091123.tsv", sep = "\t") # last version
# gr_gene_database_raw <- read.csv("data/overlapping-gr-genes/gr_geneBase_v11_141223.tsv", sep = "\t") 
gr_gene_database_raw <- read.csv("data/overlapping-gr-genes/gr_geneBase_v12_191223.tsv", sep = "\t") 
gr_gene_database_raw <- read.csv("data/overlapping-gr-genes/gr_geneBase_v15_280224.tsv", sep = "\t")
gr_gene_database_raw <- read.csv("data/overlapping-gr-genes/gr_geneBase_v16_120424.tsv", sep = "\t")

################################################################################
# Extract specific metadata from the "Info" column
gr_gene_database <- gr_gene_database_raw %>%  
  extract_keys_values("info", c("tissue", "cell", "species", "environment", "treatment", "dose", "time", "log2ratio", "fdr", "statistical_method", "treatment_type", "regulation", "comparison"))

# Filter data based on specific conditions for "source", "statistical_method", and "fdr"
tmp_michkor <- gr_gene_database %>% 
  filter(
    source == "michkor-cells", 
    statistical_method == "ttest", 
    fdr < 0.05
  )

# Merge the filtered data with the main database, excluding rows from the original "michkor-cells" source
gr_gene_database <- gr_gene_database %>% 
  filter(source != "michkor-cells") %>% 
  bind_rows(tmp_michkor)


################################################################################
# prepare hgnc_symbols from biomart v110

# Specify the URL for Ensembl version 110 from the archive
# in the future it'd change to archive url
ensembl_url_v109 = "http://feb2023.archive.ensembl.org/"
ensembl_url_v110 = "https://jul2023.archive.ensembl.org/"

# Use the biomaRt package to connect to the Ensembl dataset from the specified version (110)
# Set the biomart to "ensembl", dataset to human genes ("hsapiens_gene_ensembl"), and use the archive URL as the host
ensembl_mart_v110 <- biomaRt::useMart(biomart = "ensembl",
                                      dataset = "hsapiens_gene_ensembl",
                                      host = ensembl_url_v110)

# Query the Ensembl database for protein-coding genes and their corresponding HGNC symbols
# - attributes: The attributes we want to retrieve, i.e., the HGNC symbol and gene biotype
# - filters: We want to filter the results based on gene biotype
# - values: Specifically, we're interested in genes with a biotype of "protein_coding"
# - mart: The specific Ensembl version (110 in this case) to use for the query
hgnc_symbols_df_v110 <- biomaRt::getBM(attributes = c("hgnc_symbol"),
                                       filters = "biotype",
                                       values = "protein_coding",
                                       mart = ensembl_mart_v110)

hgnc_symbols_vector_v110 <- hgnc_symbols_df_v110$hgnc_symbol

# List available attributes for the dataset
available_attributes <- biomaRt::listAttributes(ensembl_mart_v110)
available_attributes %>% 
  filter(grepl("hgnc", name))

# preprocess the gr_gene_database
gr_gene_database %>% 
  filter(!is.na(hgnc_symbol)) %>% # remove NA hgnc_symbol
  filter(hgnc_symbol %in% hgnc_symbols_df_v110$hgnc_symbol) %>% # filter protein_coding hgnc_sybol
  # filter((hgnc_symbol %in% hgnc_to_remove)) %>% 
  mutate(gene_list_index = gsub("marpiech_.*", "marpiech_tissues_dex", gene_list_index)) %>% 
  mutate(source = ifelse(gene_list_index == "marpiech_tissues_dex", "marpiech_tissues_dex", source)) -> gr_gene_database_preproccesing



# tmp code to repair data for PMID:24777604
gr_gene_database_preproccesing %>% 
  filter(source == "pmid:24777604") %>% 
  mutate(comparison = ifelse(comparison == "Cav1-KO-DEX_vs_C57-Ethanol", "Cav1-KO-Dex_vs_Cav1-KO-Ethanol",
                             ifelse(comparison == "C57-DEX_vs_C57-Ethanol", "C57-Dex_vs_C57-Ethanol", "Cav1-KO-Ethanol_vs_C57BL6-Dex"))) %>% 
  mutate(cell = ifelse(comparison == "C57-Dex_vs_C57-Ethanol", "NPSCs", 
                       ifelse(comparison == "Cav1-KO-Dex_vs_Cav1-KO-Ethanol", "NPSCs-KO-Cav1", "NPSCs|NPSCs-KO-Cav1")),
         treatment = ifelse(comparison == "C57-Dex_vs_C57-Ethanol", "dexamethasone", 
                            ifelse(comparison == "Cav1-KO-Dex_vs_Cav1-KO-Ethanol", "dexamethasone", "vehicle-ethanol"))) -> tmp_npscs_kocav1


gr_gene_database_preproccesing %>%
  filter(source != "pmid:24777604") %>%
  rbind(., tmp_npscs_kocav1) -> gr_gene_database_preproccesing

gr_gene_database_preproccesing %>%
  filter(!grepl("omicspred_metabolon", source)) %>%
  filter(!grepl("omicspred_nithingale_", source)) %>%
  mutate(label = ifelse(
    source == "marpiech_tissues_dex",
    paste0(source, "_", gene_list_number),
    paste0(source, "_", tissue, "_", cell)
  )) -> papers_data_preprocessing


papers_data_preprocessing %>% 
  categorize_tissue(tissue_column = "tissue") %>% 
  mutate(source = ifelse(source == "pmid:22673229e", "pmid:22673229", source)) %>% 
  mutate(source = ifelse(source == "pmid:NA", "marpiech_tissues", source)) %>%  
  mutate(dose = str_remove(dose, " x")) %>% 
  mutate(dose = str_remove(dose, " 18h| 6h| 4h| 12h| 24h| 2h| 1h")) -> papers_data_preprocessing


################################################################################
# filter hgnc symbols to remove
gr_gene_database %>% 
  filter(grepl("marpiech", gene_list_index)) %>% 
  filter(gene_list_number != 884) %>% 
  # filter(source != "aall_significant_genes_marpiech") %>% 
  select(-index) %>%  
  unique() %>% 
  select(c(gene_list_number, gene_name, hgnc_symbol)) %>% 
  filter(!is.na(hgnc_symbol)) %>% 
  filter(!(gene_list_number %in% c(17,18))) %>% 
  select(-gene_name) %>% unique() %>% 
  .$hgnc_symbol %>% 
  table %>% 
  as.data.frame() %>% 
  arrange(Freq) %>%
  set_colnames(c("hgnc_symbol", "freq")) %>% 
  filter(freq > 1) %>% 
  as.data.frame() %>% .$hgnc_symbol %>% as.character() -> hgnc_to_remove

papers_data_preprocessing %>% 
  filter((source %in% c("marpiech_tissues_dex",  "marpiech_tissues"))) %>% 
  filter(!(hgnc_symbol %in% hgnc_to_remove)) %>% 
  mutate(source = ifelse(source == "marpiech_tissues_dex", "marpiech_clusters_dex", "marpiech_tissues_dex")) -> marpiech_data_preprocessing


################################################################################

  


papers_data_preprocessing %>% 
  # filter(!(treatment %in% c("TNF", "LPS", "vehicle-ethanol", "TNFalpha"))) %>% 
  filter(!(comparison %in% c( "FS30_vs_BLAM", "FS120_vs_BLAM", "FS360_vs_BLAM", "FS180_vs_BLAM"))) %>%
  # filter(dose != "0mg/kg") %>% 
  filter(source != "marpiech_tissues") %>% 
  extract_keys_values(., "info", keys = "method") %>% 
  mutate(simple_tissue = ifelse(source == "michkor-cells", "brain", simple_tissue)) %>% 
  mutate(regulation = ifelse(regulation == "dwon", "down", regulation)) %>% 
  mutate(regulation = ifelse(source == "pmid:23303060" & log2ratio > 0, "up", 
                             ifelse(source == "pmid:23303060" & log2ratio <= 0, "down", regulation))) %>% 
  filter(!(source == "michkor-cells" & hgnc_symbol == "ARHGAP8")) %>% 
  filter(!(source %in% c("marpiech_tissues_dex"))) -> papers_data_preprocessing


papers_data_preprocessing %>% 
  filter(!(source %in% c("marpiech_tissues_dex"))) %>%
  select(c(source, tissue, cell, system, simple_tissue)) %>% 
  unique %>% as.data.frame() %>% 
  write.table(file = "results/table-source-tissue-cell.tsv", sep = "\t", row.names = F, col.names = T, quote = F)



# split(papers_data_preprocessing$hgnc_symbol, papers_data_preprocessing$label) %>% lapply(., unique) -> papers_gene_list