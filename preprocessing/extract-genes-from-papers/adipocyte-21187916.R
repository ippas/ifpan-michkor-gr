# Load necessary packages and functions
source("preprocessing/preprocessing-genes-papers/install-load-packages.R")

# Read sheet names and load each sheet into a list of data frames
sheets <- excel_sheets("data/supplement-genes-papers/adipocytes-21187916/pone.0015188.s002.xls")
data <- lapply(sheets, function(sheet) read_excel("data/supplement-genes-papers/adipocytes-21187916/pone.0015188.s002.xls", sheet = sheet, skip = 1))


data[[1]] %>% 
  as.data.frame() %>% 
  .[, c(4,8:11)]  %>% 
  set_colnames(c("gene_name", "ETOH", "Dex", "FC", "fdr")) %>% 
  mutate(log2ratio = log2(FC)) %>% 
  mutate(regulation = ifelse(log2ratio > 0, "up", "down")) %>% 
  select(-c(ETOH, FC)) -> preprocessing_df

data.frame(
  article_source = "https://pubmed.ncbi.nlm.nih.gov/21187916/",
  pmid = 21187916,
  gene_name = preprocessing_df$gene_name,
  ensembl_id = NA,  # Ensembl IDs are not available
  log2ratio = preprocessing_df$log2ratio,
  regulation = preprocessing_df$regulation,
  species = "mouse",
  tissue = "adipocyte-tissue",
  cell = "3T3-L1",
  environment = "in-vitro",
  treatment = "dexamethasone",
  dose = "500nM",
  time = "6h",
  fdr = preprocessing_df$fdr,
  fdr_threshold = 0.05,  # The threshold for False Discovery Rate
  method = "microarray",   # The method used for gene expression analysis
  comparison = "dex_vs_vehicle-etanol",  # The conditions being compared
  treatment_type = "acute"  # The type of treatment
) %>% 
  write.table(
    .,
    file = "data/supplement-genes-papers/preprocessing-genes-papers/adipocyte-21187916.tsv",
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )

