# Load necessary packages and functions
source("preprocessing/preprocessing-genes-papers/install-load-packages.R")

# Load the first three sheets from the Excel file into a list of data frames, skipping the first 3 rows
sheets <- excel_sheets("data/supplement-genes-papers/prefrontal-cortex/12035_2017_586_MOESM2_ESM.xlsx")[1:3]
data <- lapply(sheets, function(sheet) read_excel("data/supplement-genes-papers/prefrontal-cortex/12035_2017_586_MOESM2_ESM.xlsx", sheet = sheet, skip = 3))

# Rename data frames in the list for clarity
names(data) <- c("S10_vs_Control", "S30_vs_Control", "S30_vs_S10")

# Process data frames to extract relevant columns, rename them, and add 'comparison' field
prefrontal_cortex_df <- imap_dfr(data, function(df, comparison_name) {
  df %>% 
    select(-baseMean) %>%  # Remove the 'baseMean' column
    set_names(c("gene_name", "ensembl_id", "log2ratio", "pvalue", "fdr")) %>%
    mutate(comparison = comparison_name)
}) %>% 
  as.data.frame() %>%
  mutate(
    # Add 'time' based on 'comparison'
    time = case_when(
      comparison == "S10_vs_Control" ~ "240h",
      comparison == "S30_vs_Control" ~ "720h",
      TRUE ~ "240h_vs_720h"
    ),
    # Add 'regulation' based on 'log2ratio'
    regulation = ifelse(log2ratio > 0, "up", "down")
  )

# Create the final database with relevant metadata and results
prefrontal_cortex_database <- data.frame(
  article_source = "https://pubmed.ncbi.nlm.nih.gov/28500512/",
  pmid = 28500512,
  gene_name = prefrontal_cortex_df$gene_name,
  ensembl_id = prefrontal_cortex_df$ensembl_id,
  log2ratio = prefrontal_cortex_df$log2ratio,
  regulation = prefrontal_cortex_df$regulation,
  species = "mouse",
  tissue = "prefrontal-cortex",
  cell = NA,
  environment = "in-vivo",
  treatment = "social-defeat-stress",
  dose = NA,
  time = prefrontal_cortex_df$time,
  pvalue = prefrontal_cortex_df$pvalue,
  fdr = prefrontal_cortex_df$fdr,
  fdr_threshold = 0.05,
  method = "RNA-seq",
  comparison = prefrontal_cortex_df$comparison,
  strain = "C57BL/6|ICR",
  statistical_method = "DESeq2",
  treatment_type = "stress_induction"
)
