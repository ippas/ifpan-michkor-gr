# Load necessary packages and functions
source("preprocessing/preprocessing-genes-papers/install-load-packages.R")

# Read sheet names and load each sheet into a list of data frames
sheets <- excel_sheets("data/supplement-genes-papers/kidney-23593176/pone.0060213.s010.xlsx")
data <- lapply(sheets, function(sheet) read_excel("data/supplement-genes-papers/kidney-23593176/pone.0060213.s010.xlsx", sheet = sheet))

names(data) <- sheets

data[[1]] %>% 
  head %>% 
  as.data.frame() %>% 
  mutate(log2FC = log(`Fold Change`))


lapply(names(data), function(name) {
  df <- data[[name]]
  df$comparison <- name
  return(df)
}) %>% 
  bind_rows() %>% 
  as.data.frame() %>% 
  filter(comparison != "Cluster of (i)") %>% 
  mutate(log2ratio = log2(`Fold Change`)) %>% 
  select(Symbol, adj.P.Val, comparison, log2ratio) %>%  
  set_colnames(c("gene_name", "fdr", "comparison", "log2ratio")) %>% 
  mutate(regulation = ifelse(log2ratio > 0, "up", "down")) %>% 
  mutate(time = "72h") %>% 
  mutate(comparison2 = comparison) %>% 
  separate(col = comparison2, 
           into = c("cell1", "treatment1", "cell2", "treatment2"), 
           sep = "-|\\.") %>% 
  # select(comparison) %>% unique() 
  mutate(cell1 = ifelse(cell1 == "D", "HPCs-differentiated", "HPCs-undifferentiated")) %>% 
  mutate(cell2 = ifelse(cell2 == "D", "HPCs-differentiated", "HPCs-undifferentiated")) %>% 
  mutate(treatment1 = ifelse(treatment1 == "Dex", "dexamethasone", ifelse(treatment1 == "VD3", "vitamin-d3", "vehicle-DMSO"))) %>% 
  mutate(treatment2 = ifelse(treatment2 == "Dex", "dexamethasone", ifelse(treatment2 == "VD3", "vitamin-d3", "vehicle-DMSO"))) %>% 
  mutate(treatment = ifelse(treatment1 == treatment2, treatment1, paste0(treatment1, "|", treatment2))) %>% 
  mutate(cell = ifelse(cell1 == cell2, cell1, paste0(cell1, "|", cell2))) %>% 
  select(-c(cell1, cell2, treatment1, treatment2)) %>%
  mutate(dose = ifelse(treatment == "vehicle-DMSO", NA, "100nM")) %>% 
  mutate(treatment_type = ifelse(treatment == "vehicle-DMSO", NA, "acute")) -> preprocessing_df


data.frame(
  article_source = "https://pubmed.ncbi.nlm.nih.gov/23593176/",
  pmid = 23593176,
  gene_name = preprocessing_df$gene_name,
  ensembl_id = NA,  # Ensembl IDs are not available
  log2ratio = preprocessing_df$log2ratio,
  regulation = preprocessing_df$regulation,
  species = "human",
  tissue = "kidney",
  cell = preprocessing_df$cell,
  environment = "in-vitro",
  treatment = preprocessing_df$treatment,
  dose = preprocessing_df$dose,
  time = preprocessing_df$time,
  fdr = preprocessing_df$fdr,
  fdr_threshold = 0.05,  # The threshold for False Discovery Rate
  method = "microarray",   # The method used for gene expression analysis
  comparison = preprocessing_df$comparison,  # The conditions being compared
  treatment_type = preprocessing_df$treatment_type  # The type of treatment
) %>% 
  write.table(
    .,
    file = "data/supplement-genes-papers/preprocessing-genes-papers//kidney-23593176.tsv",
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )



  
  


