# Load required packages and functions
source("preprocessing/preprocessing-genes-papers/install-load-packages.R")

# Define the file path of the Excel file to be read
file_path <- "data/supplement-genes-papers/placenta-PNSS/placenta-raw.xlsx"

# Retrieve the names of all sheets in the Excel file
sheet_names <- excel_sheets(file_path)

# Read the first sheet to initialize combined_data and to fetch column names
combined_data <- read_excel(file_path, sheet = sheet_names[1], skip = 1, range = cell_cols("B:Z"))
col_names <- colnames(combined_data)

# Loop over the rest of the sheets to append them to combined_data
for (i in 2:length(sheet_names)) {
  # Read the current sheet into a temporary data frame
  temp_data <- read_excel(file_path, sheet = sheet_names[i], skip = 1, range = cell_cols("B:Z"))
  
  # Make sure the column names of temp_data align with those of combined_data
  colnames(temp_data) <- col_names
  
  # Append temp_data to combined_data
  combined_data <- rbind(combined_data, temp_data)
}

# Write the combined data to a new Excel file
write_xlsx(combined_data, "data/supplement-genes-papers/placenta-PNSS/placenta-combined-data.xlsx")

# Preprocess the combined data
combined_data %>% 
  .[, 1:7] %>% # Select first 7 columns
  set_colnames(c("gene_name", "gene_type", "cell_type", "fpkm", "log2FC", "pval", "fdr")) %>% # Rename columns
  as.data.frame() %>% # Convert to data frame
  mutate(log2FC = as.numeric(log2FC)) %>% # Convert log2FC to numeric
  filter(fdr < 0.05) %>% # Filter rows by FDR < 0.05
  filter(log2FC > 1 | log2FC < -1) %>% # Filter rows by absolute log2FC > 1
  mutate(regulation = ifelse(log2FC > 0, "up", "down")) %>% # Add regulation column
  select(-c(fpkm, cell_type)) -> placenta_preprocessing # Remove fpkm and cell_type columns

# Create final data frame for the database
placanta_to_database <- data.frame(
  article_source = "https://www.nature.com/articles/s41380-021-01123-z#Sec16",
  pmid =  33981007,
  gene_name = placenta_preprocessing$gene_name,
  ensembl_id = NA,
  log2ratio = placenta_preprocessing$log2FC,
  regulation = placenta_preprocessing$regulation,
  species = "human",
  tissue = "placenta",
  cell = "CTB|STR|EVT|STB",
  environment = "in-vivo",
  treatment = "PNSS",
  dose = NA,
  time = NA, 
  pvalue = placenta_preprocessing$pval,
  fdr = placenta_preprocessing$fdr,
  fdr_threshold = 0.05,
  method = "RNA-seq",
  genome_element = placenta_preprocessing$gene_type,
  comparison = "PNSS-Negative_vs_PNSS-Positive",
  statistical_method = "DESeq2",
  treatment_type = "stress_induction"
)
