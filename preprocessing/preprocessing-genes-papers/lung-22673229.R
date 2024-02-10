source("preprocessing/functions/R/install-load-packages.R")

# read data
data_path <- "data/supplement-genes-papers/lung-22673229/supp_en.2012-1020_EN-12-1020Supp_Table_1-EL-v1.xlsx"

read_excel_sheets(file_path = data_path, skip_rows = 1) -> data

data %>% lapply(., dim)


data[[1]] %>% 
  set_colnames(c("gene_name", "log2ratio", "fdr")) %>% 
  mutate(log2ratio = as.numeric(log2ratio)) %>% 
  mutate(fdr = as.numeric(fdr)) %>% 
  na.omit() %>% 
  mutate(treatment = "dexamethasone") %>% 
  mutate(dose = "10nM") -> df_preprocessing1

data[[2]] %>%
  .[, c(1,2,3)] %>% 
  set_colnames(c("gene_name", "log2ratio", "fdr")) %>% 
  mutate(log2ratio = as.numeric(log2ratio)) %>% 
  mutate(fdr = as.numeric(fdr)) %>% 
  na.omit() %>% 
  mutate(treatment = "TNFalpha") %>% 
  mutate(dose = "10μg/ml") -> df_preprocessing2

rbind(df_preprocessing1, df_preprocessing2) %>% mutate(regulation = ifelse(log2ratio > 0, "up", "down")) -> df_preprocessing

data.frame(
  article_source = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3404340/",
  pmid = "22673229e",
  gene_name = df_preprocessing$gene_name,  # Replace with actual column name
  ensembl_id = NA,  # Replace NA with actual data if available
  log2ratio = df_preprocessing$log2ratio,  # Replace with actual column name
  regulation = df_preprocessing$regulation,  # Replace with actual column name
  species = "human",
  tissue = "lung-inflammatory",
  cell = "A549",
  environment = "in-vitro",
  treatment = df_preprocessing$treatment,
  dose = df_preprocessing$dose,  # Replace NA with actual data if available
  time = NA,  # Replace NA with actual data if available
  fdr = df_preprocessing$fdr,  # Replace with actual column name
  fdr_threshold = "0.1",
  abs_log2ratio_threshold = NA,
  method = "microarray",
  statistical_method = "ANOVA",
  treatment_type = "acute"
) %>% 
  write.table(
    .,
    file = "data/supplement-genes-papers/preprocessing-genes-papers/lung-22673229.tsv",
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )

