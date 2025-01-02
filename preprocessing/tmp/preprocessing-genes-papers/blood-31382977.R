# Load necessary packages and functions
source("preprocessing/preprocessing-genes-papers/install-load-packages.R")

data <- read.csv("data/supplement-genes-papers/blood-31382977/table-s3.csv", sep = ",", header = T)

data %>% 
  select(-c(Gene.Name, P.value)) %>% 
  set_colnames(c("gene_name", "fold_change", "fdr", "regulation")) -> preprocessing_df


data.frame(
  article_source = "https://pubmed.ncbi.nlm.nih.gov/31382977/",
  pmid = 31382977,
  gene_name = preprocessing_df$gene_name,
  ensembl_id = NA,  # Ensembl IDs are not available
  log2ratio = NA,
  regulation = preprocessing_df$regulation,
  species = "human",
  tissue = "blood-COPD",
  cell = NA,
  environment = "in-vivo",
  treatment = "prednisone",
  dose = "30mg/day",
  time = "96h",
  fdr = preprocessing_df$fdr,
  fdr_threshold = 0.05,  # The threshold for False Discovery Rate
  method = "microarray",   # The method used for gene expression analysis
  statistical_method = "t-test",
  treatment_type = "acute"  # The type of treatment
) %>%
  write.table(
    .,
    file = "data/supplement-genes-papers/preprocessing-genes-papers/blood-31382977.tsv",
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  ) 
