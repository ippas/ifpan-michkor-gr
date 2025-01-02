# Load necessary packages and functions
source("preprocessing/preprocessing-genes-papers/install-load-packages.R")

data <- read.csv("data/supplement-genes-papers/chondrocytes-33203447/13075_2020_2289_MOESM2_ESM.tsv", sep = "\t", header = T)

data %>% 
  filter(FDR < 0.05) %>% 
  mutate(abs_fc = abs(Fold_change)) %>% 
  filter(abs_fc > 2) %>% 
  mutate(log2ratio = log2(abs_fc)) %>% 
  mutate(log2ratio = ifelse(Fold_change > 0, log2ratio, (-1)*log2ratio)) %>% 
  mutate(regulation = ifelse(log2ratio > 0, "up", "down")) %>% 
  select(c(Gene, log2ratio, FDR, regulation)) %>% 
  set_colnames(c("gene_name", "log2ratiom", "fdr", "regulation")) -> preprocessing_df

data.frame(
  article_source = "https://pubmed.ncbi.nlm.nih.gov/33203447/",
  pmid = 33203447,
  gene_name = preprocessing_df$gene_name,
  ensembl_id = NA,  # Ensembl IDs are not available
  log2ratio = preprocessing_df$log2ratio,
  regulation = preprocessing_df$regulation,
  species = "human",
  tissue = "knee-cartilage",
  cell = "osteoarthritis-chon",
  environment = "in-vitro",
  treatment = "dexamethasone",
  dose = "1µM",
  time = "24h",
  fdr = preprocessing_df$fdr,
  fdr_threshold = 0.05,  # The threshold for False Discovery Rate
  method = "RNA-seq",   # The method used for gene expression analysis
  statistical_method = "DESeq2",
  treatment_type = "acute"  # The type of treatment
) %>%
  write.table(
    .,
    file = "data/supplement-genes-papers/preprocessing-genes-papers/chondrocytes-33203447.tsv",
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
)
