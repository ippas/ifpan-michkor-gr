source("preprocessing/functions/R/install-load-packages.R")

data_path <- "data/supplement-genes-papers/macrophages-37217509/41467_2023_38456_MOESM5_ESM.xlsx"

read_excel_sheets(file_path = data_path, skip_rows = 1) -> data

names(data) <- c("REH-overexpression-GCR", "NALM6")

data %>% 
  bind_rows(., .id = "cell") %>% 
  select(-baseMean) %>% 
  set_colnames(c("cell", "ensembl_id", "gene_name", "log2ratio", "pvalue", "fdr")) %>% 
  mutate(regulation = ifelse(log2ratio > 0, "up", "down")) -> df_preprocessing

data.frame(
  article_source = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10203345/",
  pmid = "37217509",
  gene_name = df_preprocessing$gene_name,  # Replace with actual column name
  ensembl_id = df_preprocessing$ensembl_id,  # Replace NA with actual data if available
  log2ratio = df_preprocessing$log2ratio,  # Replace with actual column name
  regulation = df_preprocessing$regulation,  # Replace with actual column name
  species = "human",
  tissue = "blood",
  cell = df_preprocessing$cell,
  environment = "in-vitro",
  treatment = "dexamethasone",
  dose = "1µM",  # Replace NA with actual data if available
  time = "48h",  # Replace NA with actual data if available
  pvalue = df_preprocessing$pvalue,  # Replace with actual column name
  fdr = df_preprocessing$fdr,  # Replace with actual column name
  fdr_threshold = "0.05",
  abs_log2ratio_threshold = 1,
  method = "RNA-seq",
  statistical_method = "DESeq2",
  treatment_type = "acute",
  geo = "GSE214319"
)  %>% 
  write.table(
    .,
    file = "data/supplement-genes-papers/preprocessing-genes-papers/macrophages-37217509.tsv",
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )


