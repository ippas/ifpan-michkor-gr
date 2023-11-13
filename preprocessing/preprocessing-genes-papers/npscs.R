# Load required packages and functions
source("preprocessing/preprocessing-genes-papers/install-load-packages.R")

# Read the TSV file and perform initial data transformation
# Adding a 'regulation' column based on the value of 'log2ratio'
read.table("data/supplement-genes-papers/NPSCs/npscs.tsv", header = TRUE) %>% 
  mutate(regulation = ifelse(log2ratio > 0, "up", "down")) -> npscs_df

# Create the final data frame for the database
# This includes metadata such as article source, PMIDs, and experimental conditions
data.frame(
  article_source = "https://academic.oup.com/mend/article/30/1/144/2526392",
  pmid = 26606517,
  gene_name = npscs_df$gene_name,
  ensembl_id = NA,  # Ensembl IDs are not available
  log2ratio = npscs_df$log2ratio,
  regulation = npscs_df$regulation,
  species = "mouse",
  tissue = "embryos_hypothalamic-region",
  cell = "NPSCs",
  environment = "in-vitro",
  treatment = "dexamethasone", 
  dose = "100nM",
  time = "4h",
  fdr_threshold = 0.05,  # The threshold for False Discovery Rate
  method = "RNA-seq",  # The method used for gene expression analysis
  comparison = "DEX_vs_Ethanol",  # The conditions being compared
  strain = "C57BL/6",  # The mouse strain used
  treatment_type = "acute"  # The type of treatment
) -> npscs_database  # Store the final data frame in 'npscs_database'
