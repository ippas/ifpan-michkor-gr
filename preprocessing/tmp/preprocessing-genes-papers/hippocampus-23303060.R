source("preprocessing/functions/R/install-load-packages.R")

read.csv("data/supplement-genes-papers/hippocampus-23303060/npp2012253x10.csv") %>%
  select(-Gene.Title) %>% 
  set_colnames(c("gene_name", "pvalue", "log2ratio", "dose")) %>% 
  mutate(regulation = ifelse(log2ratio > log2ratio, "up", "down")) -> df_preprocessing


data.frame(
  article_source = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3672002/",
  pmid = "23303060",
  gene_name = df_preprocessing$gene_name,  # Replace with actual column name
  ensembl_id = NA,  # Replace NA with actual data if available
  log2ratio = df_preprocessing$log2ratio,  # Replace with actual column name
  regulation = df_preprocessing$regulation,  # Replace with actual column name
  species = "human",
  tissue = "hippocampus",
  cell = "HPCO3A/07",
  environment = "in-vitro",
  treatment = "cortisol",
  dose = df_preprocessing$dose,  # Replace NA with actual data if available
  time = "10days",  # Replace NA with actual data if available
  pvalue = df_preprocessing$pvalue,  # Replace with actual column name
  fdr = NA,  # Replace with actual column name
  fdr_threshold = NA,
  abs_log2ratio_threshold = 1.2,
  method = "microarray",
  statistical_method = "ANOVA",
  treatment_type = "chronic"
) %>% 
  write.table(
    .,
    file = "data/supplement-genes-papers/preprocessing-genes-papers/hippocampus-23303060.tsv",
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )


read.csv("data/supplement-genes-papers/hippocampus-23303060/rat-hippocampus-pns.csv") %>% 
  mutate(regulation = ifelse(Group == "A", "up", "down")) %>% 
  select(-Group) %>% 
  set_colnames(c("gene_name", "pvalue", "regulation"))-> df_preprocessing

data.frame(
  article_source = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3672002/",
  pmid = "23303060",
  gene_name = df_preprocessing$gene_name,  # Replace with actual column name
  ensembl_id = NA,  # Replace NA with actual data if available
  log2ratio = NA,  # Replace with actual column name
  regulation = df_preprocessing$regulation,  # Replace with actual column name
  species = "rat",
  tissue = "hippocampus",
  cell = NA,
  environment = "in-vivo",
  treatment = "PNS",
  dose = NA,  # Replace NA with actual data if available
  time = NA,  # Replace NA with actual data if available
  pvalue = df_preprocessing$pvalue,  # Replace with actual column name
  fdr = NA,  # Replace with actual column name
  fdr_threshold = NA,
  abs_log2ratio_threshold = 1.2,
  method = "microarray",
  statistical_method = NA,
  treatment_type = "chronic"
) %>% 
  write.table(
    .,
    file = "data/supplement-genes-papers/preprocessing-genes-papers/hippocampus-rat-23303060.tsv",
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )

