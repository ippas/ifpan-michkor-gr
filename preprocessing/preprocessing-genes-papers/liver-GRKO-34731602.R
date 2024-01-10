# Load necessary packages and functions
source("preprocessing/preprocessing-genes-papers/install-load-packages.R")

read.csv("data/supplement-genes-papers/liver-GRKO-34731602/mmc2.csv") %>% 
  select(c(Gene, cluster_new)) %>% 
  set_colnames(c("gene_name", "regulation")) -> preprocessing_df


data.frame(
  article_source = "https://pubmed.ncbi.nlm.nih.gov/34731602/",
  pmid = 34731602,
  gene_name = preprocessing_df$gene_name,
  ensembl_id = NA,
  regulation = preprocessing_df$regulation,
  species = "mouse",
  tissue = "liver",
  cell = NA,
  environment = "in-vivo",
  treatment = "GRKO",
  dose = NA,
  time = "240-528h",
  fdr_threshold = 0.05,  # The threshold for False Discovery Rate
  method = "RNA-seq"   # The method used for gene expression analysis
  # comparison = "dex_vs_vehicle-etanol",  # The conditions being compared
  # treatment_type = "acute"  # The type of treatment
)  %>% 
  write.table(
    .,
    file = "data/supplement-genes-papers/preprocessing-genes-papers/liver-GRKO-34731602.tsv",
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )
