source("preprocessing/preprocessing-genes-papers/install-load-packages.R")


data_intestines <- read.csv("data/supplement-genes-papers/liver-intestine-36699046/DataSheet_1.csv", sep = ",", header = T, skip = 1)

data_intestines %>% 
  mutate(mean_veh = rowMeans(select(., Veh1, Veh2, Veh3), na.rm = TRUE)) %>%
  mutate(mean_cort = rowMeans(select(., CORT1, CORT2, CORT3), na.rm = TRUE)) %>%
  mutate(fold_change = mean_cort / mean_veh) %>% 
  mutate(log2ratio = log2(fold_change)) %>% 
  mutate(tissue = "small-intestine") %>% 
  select(-c(Veh1, Veh2, Veh3, CORT1, CORT2, CORT3, mean_veh, mean_cort, fold_change, Significance)) %>% 
  set_colnames(c("ensembl_id", "gene_name", "fdr", "regulation", "log2ratio", "tissue")) -> data_intestines
   

data_liver <- read.csv("data/supplement-genes-papers/liver-intestine-36699046/DataSheet_2.csv", sep = ",", header = T, skip = 1)

data_liver %>% 
  mutate(mean_veh = rowMeans(select(., Veh1, Veh2, Veh3), na.rm = TRUE)) %>%
  mutate(mean_cort = rowMeans(select(., CORT1, CORT2, CORT3), na.rm = TRUE)) %>%
  mutate(fold_change = mean_cort / mean_veh) %>% 
  mutate(log2ratio = log2(fold_change)) %>% 
  mutate(tissue = "liver") %>% 
  select(-c(Veh1, Veh2, Veh3, CORT1, CORT2, CORT3, mean_veh, mean_cort, fold_change, Significance)) %>% 
  set_colnames(c("ensembl_id", "gene_name", "fdr", "regulation", "log2ratio", "tissue")) -> data_liver

preprocessing_df <- rbind(data_liver, data_intestines)

data.frame(
  article_source = "https://pubmed.ncbi.nlm.nih.gov/36699046/",
  pmid = 36699046,
  gene_name = preprocessing_df$gene_name,
  ensembl_id = preprocessing_df$ensembl_id,  
  log2ratio = preprocessing_df$log2ratio,
  regulation = preprocessing_df$regulation,
  species = "mouse",
  tissue = preprocessing_df$tissue,
  cell = NA,
  environment = "in-vivo",
  treatment = "corticosterone",
  dose = "40mg/kg",
  time = "8weeks",
  fdr = preprocessing_df$fdr,
  fdr_threshold = 0.05,  # The threshold for False Discovery Rate
  method = "RNA-seq",   # The method used for gene expression analysis
  strain = "C57BL6/J",  # Strain of the animals used
  statistical_method = "DESeq2",
  treatment_type = "chronic"  # The type of treatment
) %>% 
  write.table(
    .,
    file = "data/supplement-genes-papers/preprocessing-genes-papers/liver-intestine-36699046.tsv",
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )

rm(data_liver, data_intestines)