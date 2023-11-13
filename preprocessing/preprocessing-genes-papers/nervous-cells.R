# Load necessary packages and functions
source("preprocessing/preprocessing-genes-papers/install-load-packages.R")

# Read in Excel sheets names and load each sheet into a list of data frames
sheets <- excel_sheets("data/supplement-genes-papers/nervous-cells/Gene.lists.with.annotations.xls")
data <- lapply(sheets, function(sheet) read_excel("data/supplement-genes-papers/nervous-cells/Gene.lists.with.annotations.xls", sheet = sheet))

# Standardize sheet names by converting them to lowercase, then set these as the names for the data list
sheets %>% tolower() %>% str_replace("opc", "OPC") -> sheets
names(data) <- sheets 

# Process each data frame in the list, then combine them into a single data frame
nervous_cells_df <- imap_dfr(data, ~ .x %>% 
                               # Select and rename columns
                               select(c(Symbol, Log2FoldChange, adj.P.val)) %>% 
                               set_names(c("gene_name", "log2ratio", "fdr")) %>%
                               # Add 'cell' column based on the sheet name
                               mutate(cell = .y)) %>% 
  as.data.frame() %>%
  # Filter out rows where 'gene_name' contains the word 'predicted'
  filter(!grepl("predicted", gene_name)) %>% 
  # Add 'regulation' column based on the value of 'log2ratio'
  mutate(regulation = ifelse(log2ratio > 0, "up", "down"))  # Note: Corrected "donw" to "down"

# Create the final data frame for the database
nervous_cells_database <- data.frame(
  article_source = "https://pubmed.ncbi.nlm.nih.gov/24147833/",
  pmid = 24147833,
  gene_name = nervous_cells_df$gene_name,
  ensembl_id = NA,
  log2ratio = nervous_cells_df$log2ratio,
  regulation = nervous_cells_df$regulation,
  species = "rat",
  tissue = "brain",
  cell = nervous_cells_df$cell,
  environment = "in-vitro",
  treatment = "dexamethasone",
  dose = "1μM",
  time = "48h",
  fdr = nervous_cells_df$fdr,
  fdr_threshold = 0.05,
  method = "microarray",
  treatment_type = "acute"
)
