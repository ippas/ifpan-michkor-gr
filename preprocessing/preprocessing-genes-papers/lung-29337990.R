# URL of the file
file_url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE104714&format=file&file=GSE104714%5FDEX%5FEtOH%5Fexpression%2Etsv%2Egz"

# Destination directory
dest_dir <- "data/supplement-genes-papers/lung-29337990"

# Create the directory if it doesn't exist
if (!dir.exists(dest_dir)) {
  dir.create(dest_dir, recursive = TRUE)
}

# Full path for the downloaded file
dest_file <- file.path(dest_dir, "GSE104714_DEX_EtOH_expression.tsv.gz")

# Download the file
download.file(file_url, destfile = dest_file, mode = "wb")

# Path to the downloaded file
file_path <- "data/supplement-genes-papers/lung-29337990/GSE104714_DEX_EtOH_expression.tsv.gz"

# Read the gzipped TSV file into R
expression_data <- read.table(gzfile(file_path), header = TRUE, sep = "\t",  row.names = 1)

# View the first few rows of the data
head(expression_data)

performDESeq2Analysis <- function(expression_data, pvalue_threshold = 0.05, log2fc_threshold = 2) {
  library(DESeq2)
  library(dplyr)
  
  # Create sample information data frame
  sample_info <- data.frame(
    Sample = colnames(expression_data),
    condition = c(rep("DEX", ncol(expression_data) / 2), rep("EtOH", ncol(expression_data) / 2))
  )
  
  print(sample_info)
  
  # Create DESeq2 data object
  dds <- DESeqDataSetFromMatrix(countData = expression_data,
                                colData = sample_info,
                                design = ~ condition)
  
  # Perform DESeq2 analysis
  dds <- DESeq(dds)
  
  # Get differential expression results and filter
  significant_genes <- results(dds) 
    
  
  significant_genes %>% as.data.frame %>% filter(padj < pvalue_threshold, abs(log2FoldChange) > log2(log2fc_threshold)) -> significant_genes
  
  return(significant_genes)
}

expression_data <- read.table("/home/mateusz/Downloads//GSE104714_DEX_EtOH_expression.tsv.gz", header = TRUE, sep = "\t", row.names = 1)

expression_data


time_points <- c("_1hr", "_3hr", "_5hr", "_7hr", "_9hr", "_11hr")

# Initialize an empty list to store results
results_list <- list()

for(time_point in time_points){
  filtered_expression_data <- expression_data[, grepl(time_point, colnames(expression_data))]
  significant_genes <- performDESeq2Analysis(filtered_expression_data)
  dim(significant_genes)
  time <- str_remove(time_point, "_")
  time <- str_remove(time, "r")
  results_list[[time]] <- significant_genes %>% as.data.frame() %>% rownames_to_column(var = "ensembl_id") 
  rm(significant_genes, time)
}


results_list %>% 
  imap(~ mutate(.x, time = .y)) %>%
  bind_rows() %>% 
  select(-c(baseMean, lfcSE, stat)) %>% 
  set_colnames(c("ensembl_id", "log2ratio", "pvalue", "fdr", "time")) %>% 
  mutate(regulation = ifelse(log2ratio > 0, "up", "down")) -> df_preprocessing
  

data.frame(
  article_source = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5786324/",
  pmid = "29337990",
  gene_name = NA,  # Replace with actual column name
  ensembl_id = df_preprocessing$ensembl_id,  # Replace NA with actual data if available
  log2ratio = df_preprocessing$log2ratio,  # Replace with actual column name
  regulation = df_preprocessing$regulation,  # Replace with actual column name
  species = "human",
  tissue = "lung",
  cell = "A549",
  environment = "in-vitro",
  treatment = "dexamethasone",
  dose = "100nM",  # Replace NA with actual data if available
  time = df_preprocessing$time,  # Replace NA with actual data if available
  pvalue = df_preprocessing$pvalue,  # Replace with actual column name
  fdr = df_preprocessing$fdr,  # Replace with actual column name
  fdr_threshold = 0.05,
  abs_log2ratio_threshold = 1,
  method = "RNA-seq",
  statistical_method = "DESeq2",
  treatment_type = "acute",
  geo = "GSE104714"
) %>% 
  write.table(
    .,
    file = "data/supplement-genes-papers/preprocessing-genes-papers/lung-29337990.tsv",
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )


rm(list = c("file_url", "dest_dir", "dest_file", "file_path", "expression_data", "performDESeq2Analysis"))



