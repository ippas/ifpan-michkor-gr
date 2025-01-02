# Load necessary packages and functions
source("preprocessing/preprocessing-genes-papers/install-load-packages.R")

# Load the second sheet from the first Excel file
sheets1 <- excel_sheets("data/supplement-genes-papers/hippocampus-fs/41467_2021_24967_MOESM10_ESM.xlsx")[2]
data1 <- lapply(sheets1, function(sheet) read_excel("data/supplement-genes-papers/hippocampus-fs/41467_2021_24967_MOESM10_ESM.xlsx", sheet = sheet))
names(data1) <- "exonic"

# Load the second sheet from the second Excel file
sheets2 <- excel_sheets("data/supplement-genes-papers/hippocampus-fs/41467_2021_24967_MOESM11_ESM.xlsx")[2]
data2 <- lapply(sheets2, function(sheet) read_excel("data/supplement-genes-papers/hippocampus-fs/41467_2021_24967_MOESM10_ESM.xlsx", sheet = sheet))
names(data2) <- "intronic"

# Merge and preprocess the data
hippocampus_fs_df <- imap_dfr(c(data1, data2), ~ .x %>%
                                # Rename columns
                                set_names(c("ensembl_id", "gene_name", "log2ratio", "time", "pvalue", "fdr")) %>% 
                                # Create 'comparison' variable based on 'time'
                                mutate(comparison = paste0(time, "_vs_BLAM")) %>% 
                                # Clean up 'time' values
                                mutate(time = str_replace(time, "FS", "")) %>%
                                mutate(time = str_c(time, "h")) %>% 
                                # Add 'genome_element' (exonic/intronic)
                                mutate(genome_element = .y)) %>% 
  as.data.frame() %>%
  # Add 'regulation' based on 'log2ratio'
  mutate(regulation = ifelse(log2ratio > 0, "up", "down"))

# Create the final database
hippocampus_fs_database <- data.frame(
  article_source = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8346558/",
  pmid =  34362910,
  gene_name = hippocampus_fs_df$gene_name,
  ensembl_id = hippocampus_fs_df$ensembl_id,
  log2ratio = hippocampus_fs_df$log2ratio,
  regulation = hippocampus_fs_df$regulation,
  species = "rat",
  tissue = "hippocampus",
  cell = NA,
  environment = "in-vivo",
  treatment = "force-swimming",
  dose = NA,
  time = hippocampus_fs_df$time, 
  pvalue = hippocampus_fs_df$pvalue,
  fdr = hippocampus_fs_df$fdr,
  fdr_threshold = 0.05,
  method = "RNA-seq",
  genome_element = hippocampus_fs_df$genome_element,
  comparison = hippocampus_fs_df$comparison,
  statistical_method = "edgeR",
  treatment_type = "stress_induction"
)
