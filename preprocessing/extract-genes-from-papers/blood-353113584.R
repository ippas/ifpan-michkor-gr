source("preprocessing/functions/R/install-load-packages.R")

read_tsv("data/supplement-genes-papers/blood-35313584/Differential_expression.txt", skip = 1, col_names = FALSE) %>% 
  set_colnames(c("cell", "treatment", "ensembl_id", "log2ratio", "SE", "pvalue", "fdr")) %>% 
  as.data.frame() %>% 
  mutate(abs_log2ratio = abs(log2ratio)) %>% 
  filter(abs_log2ratio > 1) %>% 
  filter(fdr < 0.1) %>% 
  mutate(regulation = ifelse(log2ratio > 0.5, "up", "down")) %>%
  mutate(dose = case_when(
    treatment == "LPS" ~ "1μg/mL",
    treatment == "PHA" ~ "2.5μg/mL",
    treatment == "LPS-DEX" ~ "1μg/mL_and_1μM",
    treatment == "PHA-DEX" ~ "2.5μg/mL_and_1μM",
    TRUE ~ treatment # Keeps the original value
  )) %>% 
  mutate(treatment = case_when(
    treatment == "LPS-DEX" ~ "LPS-dexamethasone_vs_LPS",
    treatment == "PHA-DEX" ~ "PHA-dexamethsone_vs_PHA",
    TRUE ~ treatment # Keeps the original value if none of the above conditions are met
  )) -> df_preprocessing


data.frame(
  article_source = "https://pubmed.ncbi.nlm.nih.gov/35313584/",
  pmid = "35313584",
  gene_name = NA,  # Replace with actual column name
  ensembl_id = df_preprocessing$ensembl_id,  # Replace NA with actual data if available
  log2ratio = df_preprocessing$log2ratio,  # Replace with actual column name
  regulation = df_preprocessing$regulation,  # Replace with actual column name
  species = "human",
  tissue = "blood-asthma",
  cell = df_preprocessing$cell,
  environment = "in-vitro",
  treatment = df_preprocessing$treatment,
  dose = df_preprocessing$dose,  # Replace NA with actual data if available
  time = "6h",  # Replace NA with actual data if available
  pvalue = df_preprocessing$pvalue,  # Replace with actual column name
  fdr = df_preprocessing$fdr,  # Replace with actual column name
  fdr_threshold = "0.1",
  abs_log2ratio_threshold = 0.5,
  method = "scRNA-seq",
  statistical_method = "DESeq2",
  treatment_type = "acute"
) %>% 
  write.table(
    .,
    file = "data/supplement-genes-papers/preprocessing-genes-papers/blood-35313584.tsv",
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )

