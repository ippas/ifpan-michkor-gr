# Load required packages and functions
source("preprocessing/preprocessing-genes-papers/install-load-packages.R")

# Read sheet names from the Excel file
sheets <- excel_sheets("data/supplement-genes-papers/lung/Dataset2.xls")

# Load all sheets (up to 209 rows each) from the Excel file into a list of data frames
data <- lapply(sheets, function(sheet) {
  read_excel("data/supplement-genes-papers/lung/Dataset2.xls", sheet = sheet, n_max = 209)
})

# Process the first sheet in the list (data[[1]]) to create a data frame 'lung_df'
# We're selecting specific columns and renaming them. Also, adding a 'regulation' column
data[[1]] %>% 
  select(c(gene, `Log2 Ratio`)) %>% 
  set_colnames(c("gene_name", "log2ratio")) %>% 
  as.data.frame() %>% 
  mutate(regulation = ifelse(log2ratio > 0, "up", "down")) -> lung_df

# Create the final data frame for the database, including metadata like article source, PMIDs, etc.
data.frame(
  article_source = "https://pubmed.ncbi.nlm.nih.gov/19801529/",
  pmid = 19801529,
  gene_name = lung_df$gene_name,
  ensembl_id = NA,
  log2ratio = lung_df$log2ratio,
  regulation = lung_df$regulation,
  species = "human",
  tissue = "lung",
  cell = "A549",
  environment = "in-vitro",
  treatment = "dexamethasone",
  dose = "100nM",
  time = "1h",
  fdr_threshold = 0.05,
  method = "RNA-seq",
  treatment_type = "acute"
) -> lung_database 