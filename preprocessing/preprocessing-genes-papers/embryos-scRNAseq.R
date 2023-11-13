# Load required packages and functions
source("preprocessing/preprocessing-genes-papers/install-load-packages.R")

# Load all the Excel sheets, removing some of them as specified
sheets <- excel_sheets("data/supplement-genes-papers/embryos-scRNAseq/supp_gr.276665.122_Supplemental_Table_S3.xls")[-c(1,10,11,12)]

# Read each sheet into a list of data frames
data <- lapply(sheets, function(sheet) {
  read_excel(
    "data/supplement-genes-papers/embryos-scRNAseq/supp_gr.276665.122_Supplemental_Table_S3.xls", 
    sheet = sheet, 
    skip = 1, 
    col_names = TRUE
  )
})

# Clean up sheet names for future use
sheets %>% 
  str_replace("\\.", "_") %>% 
  str_replace(" ", "") %>% 
  tolower() -> sheets

# Assign cleaned names to the data list
names(data) <- sheets 

# Separate pseudo-bulk and bulk data
data[c(3, 4, 7, 8)] -> pseudo_bulk
data[c(1,2,5,6)] -> bulk

# Process pseudo-bulk data
pseudo_bulk_df <- imap_dfr(pseudo_bulk, function(data, list_name) {
  data %>%
    select(1,2,3,6) %>%
    mutate(list_name = list_name) %>%
    set_colnames(c("gene_name", "pvalue", "fdr", "log2ratio", "list_name")) %>%  
    mutate(
      cell = str_split(list_name, "_", simplify = TRUE)[, 2],
      regulation = str_split(list_name, "_", simplify = TRUE)[, 6],
      statistical_method = "pseudo-bulk-analysis",
      fdr_threshold = 0.05
    ) %>% 
    select(-list_name)
})

# Process bulk data
bulk_df <- imap_dfr(bulk, function(data, list_name) {
  data %>%
    select(1) %>%
    mutate(
      pvalue = NA, 
      fdr = NA, 
      log2ratio = NA,
      list_name = list_name
    ) %>% 
    set_colnames(c("gene_name", "pvalue", "fdr", "log2ratio", "list_name")) %>%  
    mutate(
      cell = str_split(list_name, "_", simplify = TRUE)[, 2],
      regulation = str_split(list_name, "_", simplify = TRUE)[, 6],
      statistical_method = "bulk-analysis",
      fdr_threshold = NA
    ) %>% 
    select(-list_name)
})

# Combine pseudo-bulk and bulk data
embryos_scRNAseq_df <- rbind(pseudo_bulk_df, bulk_df)

# Create the final database with all relevant details
embryos_scRNAseq_database <- data.frame(
  article_source = "https://pubmed.ncbi.nlm.nih.gov/35948369/",
  pmid = 35948369,
  gene_name = embryos_scRNAseq_df$gene_name,
  ensembl_id = NA,
  log2ratio = embryos_scRNAseq_df$log2ratio,
  regulation = embryos_scRNAseq_df$regulation,
  species = "human", 
  tissue = "embryos",
  cell = embryos_scRNAseq_df$cell,
  environment = "in-vitro",
  treatment = "dexamethasone",
  dose = "0.1;1μg/mL",
  time = "72h",
  pvalue = embryos_scRNAseq_df$pvalue,
  fdr = embryos_scRNAseq_df$fdr,
  fdr_threshold = embryos_scRNAseq_df$fdr_threshold,
  method = "scRNAseq",
  statistical_method = embryos_scRNAseq_df$statistical_method,
  treatment_type = "chronic"
)
