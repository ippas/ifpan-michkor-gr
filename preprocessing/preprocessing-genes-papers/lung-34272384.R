source("preprocessing/functions/R/install-load-packages.R")

source("preprocessing/functions/R/preprocessing-gene-papers/download-and-read-geo-file.R")
source("preprocessing/functions/R/preprocessing-gene-papers/select_columns_by_pattern.R")
source("preprocessing/functions/R/preprocessing-gene-papers/perform_DESeq2_analysis.R")


file_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE159nnn/GSE159546/suppl/GSE159546%5Fraw%2Dcounts%2Ddexamethasone.csv.gz"
output_dir <- "data/supplement-genes-papers/lung-34272384/"

expression_data <- download_and_read_geo_file(file_url = file_url, output_dir = output_dir, separator = ";")

cells_vector <- c("A549", "H1944", "H1975", "H2122", "H460")


select_columns_by_pattern(expression_data, "A549", sort = T) %>%  
  perform_DESeq2_analysis(., condition_vector = rep(c("Hydrocortisone", "vehicle"), each = 2))  %>% 
  head

setNames(
  lapply(cells_vector, function(x) {
    select_columns_by_pattern(expression_data, pattern = x, sort_columns = TRUE) %>% 
      perform_DESeq2_analysis(., condition_vector = rep(c("Hydrocortisone", "vehicle"), each = 2))
  }), 
  cells_vector
)

lapply(cells_vector, function(x) {
  df <- select_columns_by_pattern(expression_data, pattern = x, sort_columns = TRUE) %>%
    perform_DESeq2_analysis(., condition_vector = rep(c("Hydrocortisone", "vehicle"), each = 2))
  df$cell <- x  # Add a new column with the name of the cell
  return(df)
}) %>% bind_rows() -> df_preprocessing


data.frame(
  article_source = "https://pubmed.ncbi.nlm.nih.gov/34272384/",
  pmid = "34272384",
  gene_name = NA,  # Replace with actual column name
  ensembl_id = df_preprocessing$ensembl_id,  # Replace NA with actual data if available
  log2ratio = df_preprocessing$log2ratio,  # Replace with actual column name
  regulation = df_preprocessing$regulation,  # Replace with actual column name
  species = "human",
  tissue = "lung",
  cell = df_preprocessing$cell,
  environment = "in-vitro",
  treatment = "hydrocortisone",
  dose = "2.75μM",  # Replace NA with actual data if available
  time = "8h",  # Replace NA with actual data if available
  pvalue = df_preprocessing$pvalue,  # Replace with actual column name
  fdr = df_preprocessing$fdr,  # Replace with actual column name
  fdr_threshold = "0.05",
  abs_log2ratio_threshold = 1,
  method = "RNA-seq",
  statistical_method = "DESeq2",
  treatment_type = "acute"
) %>% 
  write.table(
    .,
    file = "data/supplement-genes-papers/preprocessing-genes-papers/lung-34272384.tsv",
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )


