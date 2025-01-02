source("preprocessing/functions/R/install-load-packages.R")


file_names <- list.files(path = "data/supplement-genes-papers/lung-35902923/", pattern = "\\.txt", full.names = TRUE)[1:8]

data <- lapply(file_names, read_tsv)

names(data) <- c("ASM-siCEBPO;budesonide", "ASM-siCEBPO;TNF", "ASM-siCEBPO;budesonide_and_TNF", "ASM-siCEBPO;budesonide_and_TNF_vs_TNF",
  "ASM-siCNTR;budesonide", "ASM-siCNTR;TNF", "ASM-siCNTR;budesonide_and_TNF", "ASM-siCNTR;budesonide_and_TNF_vs_TNF")

data %>% 
  map2_df(., names())

map2_df(data, names(data), ~mutate(.x, source = .y)) %>% 
  separate(., source, into = c("cell", "treatment"), sep = ";") %>% 
  select(-c(baseMean, lfcSE, stat)) %>% 
  set_colnames(c("ensembl_id", "gene_name", "log2ratio", "pvalue", "fdr", "cell", "treatment")) %>% 
  mutate(regulation = ifelse(log2ratio > 0, "up", "down")) %>% 
  as.data.frame() %>% 
  filter(fdr < 0.05) %>% 
  mutate(abs_log2ratio = abs(log2ratio)) %>% 
  filter(abs_log2ratio > 0.58) %>%  
  mutate(dose = case_when(
    treatment == "budesonide" ~ "100nM",
    treatment == "budesonide_and_TNF" ~ "100nM_and_10ng/ml",
    treatment == "budesonide_and_TNF_vs_TNF" ~ "100nM_and_10ng/ml",
    treatment == "TNF" ~ "10ng/ml"
  )) -> df_preprocessing

data.frame(
  article_source = "https://pubmed.ncbi.nlm.nih.gov/35902923/",
  pmid = "35902923",
  gene_name = df_preprocessing$gene_name,  # Replace with actual column name
  ensembl_id = df_preprocessing$ensembl_id,  # Replace NA with actual data if available
  log2ratio = df_preprocessing$log2ratio,  # Replace with actual column name
  regulation = df_preprocessing$regulation,  # Replace with actual column name
  species = "human",
  tissue = "lung",
  cell = df_preprocessing$cell,
  environment = "in-vitro",
  treatment = df_preprocessing$treatment,
  dose = df_preprocessing$dose,  # Replace NA with actual data if available
  time = "18h",  # Replace NA with actual data if available
  pvalue = df_preprocessing$pvalue,  # Replace with actual column name
  fdr = df_preprocessing$fdr,  # Replace with actual column name
  fdr_threshold = "0.05",
  abs_log2ratio_threshold = 0.58,
  method = "RNA-seq",
  statistical_method = "DESeq2",
  treatment_type = "acute",
  geo = "GSE146017"
) %>% 
  write.table(
    .,
    file = "data/supplement-genes-papers/preprocessing-genes-papers/lung-35902923.tsv",
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )



rm(file_names)
