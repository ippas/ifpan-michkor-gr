# Load necessary packages and functions
source("preprocessing/preprocessing-genes-papers/install-load-packages.R")

read.csv("data/supplement-genes-papers/lung-33923915/Table_S1.tsv", sep = "\t",
         header = TRUE) %>% 
  select(c("gene.name", "gene.regulation", "fdr")) %>% 
  set_colnames(c("gene_name", "regulation", "fdr")) %>% 
  mutate(regulation = ifelse(regulation == "up-regulated", "up", "down")) -> preprocessing_df


data.frame(
  article_source = "https://pubmed.ncbi.nlm.nih.gov/33923915/",
  pmid = 33923915,
  gene_name = preprocessing_df$gene_name,
  ensembl_id = NA,  # Ensembl IDs are not available
  regulation = preprocessing_df$regulation,
  species = "human",
  tissue = "lung",
  cell = "A549",
  environment = "in-vitro",
  treatment = "dexamethasone",
  dose = "100nM",
  time = "0-12h",
  fdr = preprocessing_df$fdr,
  fdr_threshold = 0.0000001,  # The threshold for False Discovery Rate
  method = "RNA-seq",   # The method used for gene expression analysis
  treatment_type = "acute",
  statistical_method = "ANOVA"# The type of treatment
) %>% 
  write.table(
    .,
    file = "data/supplement-genes-papers/preprocessing-genes-papers/lung-33923915.tsv",
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )
