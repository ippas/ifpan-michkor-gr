# Load necessary packages and functions
source("preprocessing/functions/R/install-load-packages.R")

# read data
data_path <- "data/supplement-genes-papers/liver-34375333/pgen.1009737.s015.csv"

read.csv(file = data_path, header = TRUE, sep = ",") %>% 
  select(-c(baseMean, baseMean.1, lfcSE, stat, lfcSE.1, stat.1, baseMean.2, lfcSE.2, stat.2)) %>% 
  gather(key = "key", value = "value", -Ensemble) %>% 
  mutate(treatment = ifelse(key %in% c("VEH.exon.log2_FC", "pvalue", "padj"), "vehicle", "corticosterone")) %>% 
  mutate(dose = ifelse(key %in% c("VEH.exon.log2_FC", "pvalue", "padj"), "pulsative/constant", 
                       ifelse(key %in% c("PLS.exon.log2_FC", "pvalue.1", "padj.1"), "pulsative", "constant"))) %>% 
  mutate(column = case_when(
    grepl("log2_FC", key) ~ "log2_FC",
    grepl("pvalue", key) ~ "pvalue",
    grepl("padj", key) ~ "padj",
    TRUE ~ as.character(key) # Default case to keep original value if none of the above match
  )) %>% 
  select(-key) %>% 
  spread(., key = "column", value = "value") %>% 
  set_colnames(c("ensembl_id", "treatment", "dose", "log2ratio", "fdr", "pvalue")) %>% 
  na.omit() %>% 
  filter(fdr < 0.05) %>% 
  mutate(abs_log2ratio = abs(log2ratio)) %>% 
  filter(abs_log2ratio > log2(1.5)) %>%
  mutate(regulation = ifelse(log2ratio > 0, "up", "down")) %>% 
  filter(treatment != "vehicle") -> df_preprocessing


data.frame(
  article_source = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8378686/#pgen.1009737.s006",
  pmid = "34375333",
  gene_name = NA,  # Replace with actual column name
  ensembl_id = df_preprocessing$ensembl_id,  # Replace NA with actual data if available
  log2ratio = df_preprocessing$log2ratio,  # Replace with actual column name
  regulation = df_preprocessing$regulation,  # Replace with actual column name
  species = "rat",
  tissue = "liver",
  cell = NA,
  environment = "in-vivo",
  treatment = df_preprocessing$treatment,
  dose = df_preprocessing$dose,  # Replace NA with actual data if available
  time = NA,  # Replace NA with actual data if available
  pvalue = df_preprocessing$pvalue,  # Replace with actual column name
  fdr = df_preprocessing$fdr,  # Replace with actual column name
  fdr_threshold = "0.05",
  abs_log2ratio_threshold = 0.58,
  method = "RNA-seq",
  statistical_method = "DESeq2",
  treatment_type = "acute",
  geo = "GSE171647"
) %>% 
  write.table(
    .,
    file = "data/supplement-genes-papers/preprocessing-genes-papers/liver-34375333.tsv",
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )
  
