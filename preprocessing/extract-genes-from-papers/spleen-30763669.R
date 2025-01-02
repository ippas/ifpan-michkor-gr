# Load necessary packages and functions
source("preprocessing/preprocessing-genes-papers/install-load-packages.R")


# Read sheet names and load each sheet into a list of data frames
path <- "data/supplement-genes-papers/spleen-30763669/1-s2.0-S0378111919301301-mmc1.xlsx"
sheets <- excel_sheets(path)
data <- lapply(sheets, function(sheet) read_excel(path, sheet = sheet))

data[[2]] %>% 
  select(-c(CON_1, CON_2, DEX_1, DEX_2, Chrom, Strand, Start, End, GeneLength, GeneType, GeneDescription)) %>% 
  set_colnames(c("ensembl_id", "log2ratio", "pvalue", "fdr", "gene_name")) %>% 
  mutate(regulation = ifelse(log2ratio > 0, "up", "down")) -> preprocessing_df



data.frame(
  article_source = "https://pubmed.ncbi.nlm.nih.gov/30763669/",
  pmid = "3073669",
  gene_name = preprocessing_df$gene_name,  # Replace with actual column name
  ensembl_id = NA,  # Replace NA with actual data if available
  log2ratio = preprocessing_df$log2ratio,  # Replace with actual column name
  regulation = preprocessing_df$regulation,  # Replace with actual column name
  species = "mouse",
  tissue = "spleen",
  cell = NA,
  environment = "in-vivo",
  treatment = "dexamethasone",
  dose = "10µg/g/day",  # Replace NA with actual data if available
  time = "168h",  # Replace NA with actual data if available
  pvalue = preprocessing_df$pvalue,  # Replace with actual column name
  fdr = preprocessing_df$fdr,  # Replace with actual column name
  fdr_threshold = 0.05,
  method = "RNA-seq",
  strain = "C57BL/6",
  statistical_method = "DESeq",
  treatment_type = "chronic"
)  %>% 
  write.table(
    .,
    file = "data/supplement-genes-papers/preprocessing-genes-papers/spleen-30763669.tsv",
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )

