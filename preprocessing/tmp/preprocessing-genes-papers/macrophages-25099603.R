source("preprocessing/functions/R/install-load-packages.R")

data_path <- "data/supplement-genes-papers/macrophages-25099603/table-S2.tsv"

read.table(data_path, sep = "\t", header = TRUE, skip = 3) -> data
  
data %>%
  .[, c(1,2,3,4,5,6)] %>%
  na.omit() %>%
  set_colnames(c("gene_name", "fdr", "U", "D", "L", "L_and_D")) %>%
  mutate(log2ratio_D = log2(D/U)) %>%
  mutate(log2ratio_L = log2(L/U)) %>%
  mutate(log2ratio_L_and_D = log2(L_and_D/U)) %>%
  select(-c(U, D, L, L_and_D)) %>%
  gather(., key = "treatment", "log2ratio", -c(gene_name, fdr)) %>%
  filter(fdr < 0.1) %>%
  mutate(abs_log2ratio = abs(log2ratio)) %>%
  filter(abs_log2ratio > log(1.5)) %>%
  mutate(regulation = ifelse(log2ratio > 0, "up", "down")) %>%
  select(-abs_log2ratio) %>%
  mutate(
    treatment = case_when(
      treatment == "log2ratio_D" ~ "dexamethasone",
      treatment == "log2ratio_L" ~ "LPS",
      treatment == "log2ratio_L_and_D" ~ "LPS_and_dexamethasone",
      TRUE ~ treatment  # Keep the original treatment if none of the conditions are met
    )
  ) %>% 
  mutate(
    dose = case_when(
      treatment == "dexamethasone" ~ "100nM",
      treatment == "LPS" ~ "10ng/ml",
      treatment == "LPS_and_dexamethasone" ~ "10ng/ml_and_100nM",
      TRUE ~ treatment # Keep the original treatment if none of the conditions are met
    )
  ) -> df_preprocessing


data.frame(
  article_source = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4133603/",
  pmid = "25099603",
  gene_name = df_preprocessing$gene_name,  # Replace with actual column name
  ensembl_id = NA,  # Replace NA with actual data if available
  log2ratio = df_preprocessing$log2ratio,  # Replace with actual column name
  regulation = df_preprocessing$regulation,  # Replace with actual column name
  species = "mouse",
  tissue = "blood",
  cell = "macrophages",
  environment = "in-vitro",
  treatment = df_preprocessing$treatment,
  dose = df_preprocessing$dose,  # Replace NA with actual data if available
  time = "1h",  # Replace NA with actual data if available
  fdr = df_preprocessing$fdr,  # Replace with actual column name
  fdr_threshold = "0.1",
  abs_log2ratio_threshold = 0.58,
  method = "RNA-seq",
  strain = "C57BL/g",
  statistical_method = "ANOVA",
  treatment_type = "acute"
) %>% 
  write.table(
    .,
    file = "data/supplement-genes-papers/preprocessing-genes-papers/macrophages-25099603.tsv",
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )

