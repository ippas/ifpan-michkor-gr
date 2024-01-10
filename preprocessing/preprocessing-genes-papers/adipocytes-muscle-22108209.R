source("preprocessing/preprocessing-genes-papers/install-load-packages.R")

# Read sheet names and load each sheet into a list of data frames
sheets1 <- excel_sheets("data/supplement-genes-papers/adipocytes-muscle-22108209/tableS1.xls")
data1 <- lapply(sheets1, function(sheet) read_excel("data/supplement-genes-papers/adipocytes-muscle-22108209/tableS1.xls", sheet = sheet, skip = 2))


data1 %>% 
  as.data.frame() %>% 
  filter(!is.na(Gene.symbol)) %>% 
  select(-c(Gene.Name, Function)) %>% 
  set_colnames(c("gene_name", "fold_change")) %>%
  mutate(log2ratio = log2(fold_change)) %>% 
  mutate(regulation = ifelse(log2ratio > 0, "up", "down")) %>% 
  mutate(abs_log2ratio = abs(log2ratio)) %>% 
  filter(abs_log2ratio > 0.58) %>% 
  select(-c(abs_log2ratio, fold_change)) %>% 
  mutate(tissue = "adipose") -> preprocessing_df1
  
  
sheets2 <- excel_sheets("data/supplement-genes-papers/adipocytes-muscle-22108209/tableS2.xls")
data2 <- lapply(sheets2, function(sheet) read_excel("data/supplement-genes-papers/adipocytes-muscle-22108209/tableS2.xls", sheet = sheet, skip = 2))

data2 %>% 
  as.data.frame() %>% 
  filter(!is.na(Gene.Symbol)) %>% 
  select(-c(Gene.name, Function)) %>% 
  set_colnames(c("gene_name", "fold_change")) %>%
  mutate(log2ratio = log2(fold_change)) %>% 
  mutate(regulation = ifelse(log2ratio > 0, "up", "down")) %>% 
  mutate(abs_log2ratio = abs(log2ratio)) %>% 
  filter(abs_log2ratio > 0.58) %>% 
  select(-c(abs_log2ratio, fold_change)) %>% 
  mutate(tissue = "skeletal-muscle") -> preprocessing_df2

bind_rows(preprocessing_df1, preprocessing_df2) -> preprocessing_df


data.frame(
  article_source = "https://pubmed.ncbi.nlm.nih.gov/22108209/",
  pmid = 22108209,
  gene_name = preprocessing_df$gene_name,
  ensembl_id = NA,  # Ensembl IDs are not available
  log2ratio = preprocessing_df$log2ratio,
  regulation = preprocessing_df$regulation,
  species = "human",
  tissue = preprocessing_df$tissue,
  cell = NA,
  environment = "in-vivo",
  treatment = "dexamethasone",
  dose = "4mg/day",
  time = "96h",
  fdr = NA,
  fdr_threshold = 0.05,  # The threshold for False Discovery Rate
  method = "microarray",   # The method used for gene expression analysis
  treatment_type = "acute"  # The type of treatment
) %>% 
  write.table(
    .,
    file = "data/supplement-genes-papers/preprocessing-genes-papers/adipocyte-muscle-22108209.tsv",
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )

