# Load necessary packages and functions
source("preprocessing/preprocessing-genes-papers/install-load-packages.R")

# ----- First Data File -----

# Read sheet names and load each sheet into a list of data frames
sheets1 <- excel_sheets("data/supplement-genes-papers/adrenal-gland/supplementary_table_2.xlsx")
data1 <- lapply(sheets1, function(sheet) read_excel("data/supplement-genes-papers/adrenal-gland/supplementary_table_2.xlsx", sheet = sheet))

# Keep only relevant columns
data1 <- lapply(data1, function(x) x %>% select(c(gene_id, gene_name, log2FoldChange, pvalue, padj)))

# ----- Second Data File -----

# Read sheet names and load each sheet into a list of data frames
sheets2 <- excel_sheets("data/supplement-genes-papers/adrenal-gland/supplementary_table_3.xlsx")
data2 <- lapply(sheets2, function(sheet) read_excel("data/supplement-genes-papers/adrenal-gland/supplementary_table_3.xlsx", sheet = sheet))

# Rename columns and keep only relevant columns
data2 <- lapply(data2, function(x) x %>% set_colnames(c("gene_id", "gene_name", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")) %>% select(c(gene_id, gene_name, log2FoldChange, pvalue, padj)))

# Combine both sets of data
adrenal_gland_data <- c(data1, data2)

# Standardize sheet names
names(adrenal_gland_data) <- c(sheets1, sheets2) %>% gsub(" ", "_", .)

# Save the combined data to an Excel file
write_xlsx(adrenal_gland_data, "data/supplement-genes-papers/adrenal-gland/adrenal_gland.xlsx")

# Add a new column to each data frame with the name of the sheet it came from
adrenal_gland_data_named <- lapply(names(adrenal_gland_data), function(name) {
  df <- adrenal_gland_data[[name]]
  df$sheet_name <- name
  return(df)
})

# Combine all the data frames into a single one
adrenal_gland_df <- bind_rows(adrenal_gland_data_named) %>% 
  set_colnames(c("ensembl_id", "gene_name", "log2ratio", "pvalue", "fdr", "sheet_name")) %>% 
  mutate(
    # Extract 'environment', 'time', and 'regulation' based on the sheet name
    environment = str_replace_all(str_extract(sheet_name, "in_vivo|in_vitro"), "_", "-"),
    time = ifelse(environment == "in-vivo", "1H", str_extract(sheet_name, "1H|24H")),
    time = tolower(time),
    regulation = tolower(str_extract(sheet_name, "UP|DOWN|Up|Down"))
  ) %>%
  # Remove unnecessary columns and filter by FDR < 0.05
  select(-c(sheet_name)) %>% 
  filter(fdr < 0.05)

# Create the final data frame for the database
adrenal_gland_database <- data.frame(
  article_source = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9346337/",
  pmid = 35904237,
  gene_name = adrenal_gland_df$gene_name,
  ensembl_id = adrenal_gland_df$ensembl_id,
  log2ratio = adrenal_gland_df$log2ratio,
  regulation = adrenal_gland_df$regulation,
  species = "mouse",
  tissue = "adrenal-gland",
  cell = "Y-1",
  environment = adrenal_gland_df$environment,
  treatment = "dexamethasone",
  dose = "100nM",
  time = adrenal_gland_df$time,
  pvalue = adrenal_gland_df$pvalue,
  fdr = adrenal_gland_df$fdr,
  fdr_threshold = 0.05,
  method = "RNA-seq",
  statistical_method = "DESeq2",
  treatment_type = "acute"
)
