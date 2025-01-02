# Load necessary packages and functions
source("preprocessing/preprocessing-genes-papers/install-load-packages.R")

data <- read.csv("data/supplement-genes-papers/macrophages-24395918/", sep = ",", header = T)


sheets <- excel_sheets("data/supplement-genes-papers/macrophages-24395918/ji_1302138_supplemental_material_1.xlsx")
data <- lapply(sheets[[2]], function(sheet) read_excel("data/supplement-genes-papers/macrophages-24395918/ji_1302138_supplemental_material_1.xlsx", sheet = sheet))

data %>% 
  as.data.frame() %>% 
  select(c(Symbol, Fold.Change)) %>% 
  set_colnames(c("gene_name", "fold_change")) %>% 
  mutate(regulation = ifelse(fold_change > 0, "up", "down")) -> preprocessing_df


data.frame(
  article_source = "https://pubmed.ncbi.nlm.nih.gov/24395918/",
  pmid = 24395918,
  gene_name = preprocessing_df$gene_name,
  ensembl_id = NA,  # Ensembl IDs are not available
  log2ratio = NA,
  regulation = preprocessing_df$regulation,
  species = "human",
  tissue = "blood",
  cell = "macrophages",
  environment = "in-vitro",
  treatment = "fluticasone-propionate",
  dose = "100nM",
  time = "168h",
  fdr_threshold = 0.05,  # The threshold for False Discovery Rate
  method = "microarray",   # The method used for gene expression analysis
  treatment_type = "chronic"  # The type of treatment
) %>%
  write.table(
    .,
    file = "data/supplement-genes-papers/preprocessing-genes-papers/macrophages-24395918.tsv",
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  ) 
