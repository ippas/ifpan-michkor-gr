# Load necessary packages and functions
source("preprocessing/preprocessing-genes-papers/install-load-packages.R")

sheets <- excel_sheets("data/supplement-genes-papers/kidney-34735568/BSR-2021-1847_supp.xlsx")[1]
data <- lapply(sheets, function(sheet) read_excel("data/supplement-genes-papers/kidney-34735568/BSR-2021-1847_supp.xlsx", sheet = sheet))

data %>% 
  .[[1]] %>% 
  select(c(ID, log2FoldChange, padj, GeneName, GeneType)) %>% 
  set_colnames(c("ensembl_id", "log2ratio", "fdr", "gene_name", "genome_element")) %>% 
  as.data.frame() %>% 
  mutate(regulation = ifelse(log2ratio > 0, "up", "down")) %>% 
  filter(fdr < 0.05) -> preprocessing_df


  
data.frame(
  article_source = "https://pubmed.ncbi.nlm.nih.gov/34735568/",
  pmid = 34735568,
  gene_name = preprocessing_df$gene_name,
  ensembl_id = preprocessing_df$ensembl_id,  # Ensembl IDs are not available
  log2ratio = preprocessing_df$log2ratio,
  regulation = preprocessing_df$regulation,
  species = "mouse",
  tissue = "kidney",
  cell = NA,
  environment = "in-vivo",
  treatment = "dexamethasone",
  dose = "10µg/g/day",
  time = "168h",
  fdr = preprocessing_df$fdr,
  fdr_threshold = 0.05,  # The threshold for False Discovery Rate
  method = "RNA-seq",   # The method used for gene expression analysis
  strain = "C57BL/6",
  treatment_type = "chronic"  # The type of treatment
) %>%
  write.table(
    .,
    file = "data/supplement-genes-papers/preprocessing-genes-papers/kidney-34735568.tsv",
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )


