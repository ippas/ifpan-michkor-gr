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
  install.packages("biomaRt")
}

# Detach the biomaRt package
detach(package:biomaRt)

# Load in custom functions related to preprocessing
source("preprocessing/functions/R/gr-database-functions.R")
source("preprocessing/functions/R/draw-custom-heatmap.R")
source("preprocessing/functions/R/processing-overlap-results.R")
source("preprocessing/functions/R/manual-filter-overlap-results.R")
source("preprocessing/functions/R/gene-paper-preprocessing-functions.R")


# Read data
# gr_gene_database_raw <- read.csv("data/overlapping-gr-genes/geneBase_060723.tsv", sep = "\t") 
# gr_gene_database_raw <- read.csv("data/overlapping-gr-genes/publikacje_gr_100923.tsv", sep = "\t") 
# gr_gene_database_raw <- read.csv("data/overlapping-gr-genes/publikacje_gr_v5_181023.tsv", sep = "\t")
# gr_gene_database_raw <- read.csv("data/overlapping-gr-genes/gr_geneBase_v8_251023.tsv", sep = "\t") 
# gr_gene_database_raw <- read.csv("data/overlapping-gr-genes/gr_geneBase_v9_091123.tsv", sep = "\t") # last version
# gr_gene_database_raw <- read.csv("data/overlapping-gr-genes/gr_geneBase_v11_141223.tsv", sep = "\t") 
gr_gene_database_raw <- read.csv("data/overlapping-gr-genes/gr_geneBase_v12_191223.tsv", sep = "\t") 



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

# filter hgnc symbols to remove
gr_gene_database %>% 
  filter(grepl("marpiech", gene_list_index)) %>% 
  filter(gene_list_number != 800) %>% 
  filter(source != "marpiech_tissues_dex") %>% 
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
  filter(!(hgnc_symbol %in% hgnc_to_remove)) %>% 
  mutate(gene_list_index = gsub("marpiech_.*", "marpiech_tissues_dex", gene_list_index)) %>% 
  mutate(source = ifelse(gene_list_index == "marpiech_tissues_dex", "marpiech_tissues_dex", source)) -> gr_gene_database_preproccesing
# mutate(source = str_replace_all(source, "pmid:", ""))
