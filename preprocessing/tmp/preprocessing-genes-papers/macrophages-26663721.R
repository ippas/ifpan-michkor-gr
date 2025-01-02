# Load necessary packages and functions
source("preprocessing/preprocessing-genes-papers/install-load-packages.R")

# read data
data_path <- "data/supplement-genes-papers/macrophages-26663721/NIHMS65998-supplement-Supplementary_Tables.xlsx"

read_excel_sheets(file_path = data_path, skip_rows = 5) -> data

data[[1]] %>% 
  set_colnames(c("gene_name", "log2ratio", "fdr", "time")) %>% 
  mutate(regulation = ifelse(log2ratio > 0, "up", "down")) %>% 
  mutate(species = "mouse") %>% 
  mutate(strain = "C57BL/6") %>% 
  mutate(cell = "mBMDM") -> df_preprocessing1

data[[2]] %>% 
  set_colnames(c("gene_name", "log2ratio", "fdr", "time")) %>% 
  mutate(regulation = ifelse(log2ratio > 0, "up", "down")) %>% 
  mutate(species = "human") %>% 
  mutate(strain = "NA") %>% 
  mutate(cell = "hMDM") -> df_preprocessing2

rbind(df_preprocessing1, df_preprocessing2) -> df_preprocessing

rm(df_preprocessing1, df_preprocessing2)

data.frame(
  article_source = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4707550/",
  pmid = "26663721",
  gene_name = df_preprocessing$gene_name,  # Replace with actual column name
  ensembl_id = NA,  # Replace NA with actual data if available
  log2ratio = df_preprocessing$log2ratio,  # Replace with actual column name
  regulation = df_preprocessing$regulation,  # Replace with actual column name
  species = df_preprocessing$species,
  tissue = "blood",
  cell = df_preprocessing$cell,
  environment = "in-vitro",
  treatment = "dexamethasone",
  dose = "100nM",  # Replace NA with actual data if available
  time = df_preprocessing$time,  # Replace NA with actual data if available
  fdr = df_preprocessing$fdr,  # Replace with actual column name
  fdr_threshold = "0.05",
  abs_log2ratio_threshold = 1,
  method = "microarray",
  strain = df_preprocessing$strain,
  statistical_method = "arrayQualityMetrics,affy,limma",
  treatment_type = "acute",
  geo = "GSE61881"
) %>% 
  write.table(
    .,
    file = "data/supplement-genes-papers/preprocessing-genes-papers/macrophages-26663721.tsv",
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )

