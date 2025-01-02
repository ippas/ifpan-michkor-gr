source("preprocessing/functions/R/install-load-packages.R")

# read data

read_tsv("data/supplement-genes-papers/bronchial-28116096/table-S3.tsv", skip = 2) %>% 
  na.omit() %>% 
  .[, c(3,5,6)] %>% 
  set_colnames(c("gene_name", "fdr", "fold")) %>% 
  mutate(log2ratio = log2(fold)) %>% 
  mutate(tissue = "bronchial") %>% 
  mutate(cell = NA) %>% 
  mutate(environment = "in-vivo") %>% 
  mutate(dose = "1600μg") -> df_preprocessing1

read_tsv("data/supplement-genes-papers/bronchial-28116096/table-S9.tsv", skip = 2) %>% 
  na.omit() %>% 
  .[, c(3,5,6)] %>% 
  set_colnames(c("gene_name", "fdr", "fold")) %>% 
  mutate(log2ratio = log2(fold)) %>% 
  mutate(tissue = "bronchial") %>% 
  mutate(cell = "HBE") %>% 
  mutate(environment = "in-vitro") %>% 
  mutate(dose = "0.1μM") -> df_preprocessing2

rbind(df_preprocessing1, df_preprocessing2) %>% 
  mutate(regulation = ifelse(log2ratio > 0, "up", "dwon")) -> df_preprocessing


data.frame(
  article_source = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5242176/",
  pmid = "28116096",
  gene_name = df_preprocessing$gene_name,  # Replace with actual column name
  ensembl_id = NA,  # Replace NA with actual data if available
  log2ratio = df_preprocessing$log2ratio,  # Replace with actual column name
  regulation = df_preprocessing$regulation,  # Replace with actual column name
  species = "human",
  tissue = df_preprocessing$tissue,
  cell = df_preprocessing$cell,
  environment = df_preprocessing$environment,
  treatment = "budesonide",
  dose = df_preprocessing$dose,  # Replace NA with actual data if available
  time = "6h",  # Replace NA with actual data if available
  pvalue = NA,  # Replace with actual column name
  fdr = df_preprocessing$fdr,  # Replace with actual column name
  fdr_threshold = "0.05",
  abs_log2ratio_threshold = 1,
  method = "microarray",
  statistical_method = "ANOVA",
  treatment_type = "acute",
  geo = "GSE83233"
) %>% 
  write.table(
    .,
    file = "data/supplement-genes-papers/preprocessing-genes-papers/bronchial-28116096.tsv",
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )
