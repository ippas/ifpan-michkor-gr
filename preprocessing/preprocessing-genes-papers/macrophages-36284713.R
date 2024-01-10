# Load necessary packages and functions
source("preprocessing/preprocessing-genes-papers/install-load-packages.R")

# Read sheet names and load each sheet into a list of data frames
sheets <- excel_sheets("data/supplement-genes-papers/macrophage-36284713/mmc4.xlsx")
data <- lapply(sheets, function(sheet) read_excel("data/supplement-genes-papers/macrophage-36284713/mmc4.xlsx", sheet = sheet))


data[[3]] %>% 
  as.data.frame() %>% 
  .[, c(2:5)] %>% 
  set_colnames(c("ensemebl_id", "log2ratio", "fdr", "color")) %>% 
  filter(color %in% c("blue", "red")) %>% 
  mutate(regulation = ifelse(color == "red", "up", "down")) %>% 
  select(-color) ->  preprocessing_df


data.frame(
  article_source = "https://pubmed.ncbi.nlm.nih.gov/36284713/",
  pmid = 36284713,
  gene_name = NA,
  ensembl_id = preprocessing_df$ensemebl_id,  # Ensembl IDs are not available
  log2ratio = preprocessing_df$log2ratio,
  regulation = preprocessing_df$regulation,
  species = "human",
  tissue = "blood",
  cell = "THP-1",
  environment = "in-vitro",
  treatment = "dexamethasone",
  dose = "1μM",
  time = "3h",
  fdr = preprocessing_df$fdr,
  fdr_threshold = 0.05,  # The threshold for False Discovery Rate
  method = "RNA-seq",   # The method used for gene expression analysis
  comparison = "dex_vs_vehicle-etanol",  # The conditions being compared
  treatment_type = "acute",
  statistical_method = "DESeq2" # The type of treatment
) %>% 
  write.table(
    .,
    file = "data/supplement-genes-papers/preprocessing-genes-papers/macrophages-36284713.tsv",
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )
  
