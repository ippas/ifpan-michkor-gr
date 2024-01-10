# Load necessary packages and functions
source("preprocessing/preprocessing-genes-papers/install-load-packages.R")


file_path <- "data/supplement-genes-papers/skeletal-muscle-28432191/Supplementary_Table_1.tsv"

process_rows(df, column = "text, separtor")

split_columns <- function(df, column, separator) {
  max_values <- max(str_count(df[[column]], separator)) + 1
  column_names <- paste0("column", 1:max_values)
  
  split_df <- df %>%
    separate({{column}}, into = column_names, sep = separator, fill = "right")
  
  return(split_df)
}

process_rows <- function(df, column, separator) {


  asplit(df, 1) %>%  lapply(., function(x) {
    as.data.frame(x) %>%
      set_colnames("text") %>% 
      split_columns(., "text", "\n") %>%  
      t %>% 
      as.data.frame %>% 
      rownames_to_column(var = "row") %>% 
      select(-row)
  }) -> tmp_list
  
  as.data.frame(do.call(rbind, tmp_list)) -> final_df
  
  return(final_df)
}

file_path <- "data/supplement-genes-papers/skeletal-muscle-28432191/Supplementary_Table_1.tsv"

read_tsv(file_path, col_names = F) %>%
  .[-c(1:4), c(1,3,4)] %>% 
  # head %>% 
  na.omit() %>% 
  set_colnames(c("gene_name", "fc_aldosterone", "fc_prednisolone")) %>% 
  as.data.frame()  %>%
  process_rows() -> df1
  
df1 %>% 
  select(c(gene_name, fc_aldosterone)) %>% 
  filter(fc_aldosterone != "X") %>% 
  set_colnames(c("gene_name", "fc")) %>% 
  mutate(treatment = "aldosterone",
         dose = "10µM") -> df_aldosterone


df1 %>% 
  select(c(gene_name, fc_prednisolone)) %>% 
  filter(fc_prednisolone != "X") %>% 
  set_colnames(c("gene_name", "fc")) %>% 
  mutate(treatment = "prednisolone", 
         dose = "1µM")-> df_prednisolone

file_path <- "data/supplement-genes-papers/skeletal-muscle-28432191/Supplementary_Table_2.tsv"


read_tsv(file_path, col_names = F) %>%
  .[-c(1:4), c(1,3)] %>% 
  # head %>% 
  na.omit() %>% 
  set_colnames(c("gene_name", "fc_eplerenome")) %>% 
  as.data.frame()  %>% 
  process_rows() %>% 
  set_colnames(c("gene_name", "fc")) %>% 
  mutate(treatment = "eplerenone",
         dose = "10µM") -> df_eplerenone


file_path <- "data/supplement-genes-papers/skeletal-muscle-28432191/Supplementary_Table_3.tsv"
read_tsv(file_path, col_names = F) %>%
  .[-c(1:4), c(1,3)] %>% 

  na.omit() %>% 
  set_colnames(c("gene_name", "fc_mifepristone")) %>% 
  as.data.frame()  %>% 
  process_rows() %>% 
  set_colnames(c("gene_name", "fc")) %>% 
  mutate(treatment = "mifepristone", 
         dose = "1µM") -> df_mifepristone



rbind(df_prednisolone, df_aldosterone, df_eplerenone, df_mifepristone) %>% 
  mutate(fc = as.numeric(fc)) %>% 
  mutate(log2ratio = ifelse(fc > 0, log2(fc), -log2(abs(fc)))) %>% 
  mutate(regulation = ifelse(log2ratio > 0, "up", "down")) -> preprocessing_df


data.frame(
  article_source = "https://pubmed.ncbi.nlm.nih.gov/28432191/",
  pmid =  28432191,
  gene_name = preprocessing_df$gene_name,
  ensembl_id = NA,  # Ensembl IDs are not available
  log2ratio = preprocessing_df$log2ratio,
  regulation = preprocessing_df$regulation,
  species = "human",
  tissue = "skeletal-muscle",
  cell = "myotubes",
  environment = "in-vitro",
  treatment = preprocessing_df$treatment,
  dose = preprocessing_df$dose,
  time = "48h",
  fdr = NA,
  fdr_threshold = NA,  # The threshold for False Discovery Rate
  method = "microarray",   # The method used for gene expression analysis
  # comparison = "dex_vs_vehicle-etanol",  # The conditions being compared
  treatment_type = "acute"  # The type of treatment
) %>% 
  write.table(
    .,
    file = "data/supplement-genes-papers/preprocessing-genes-papers/skeletal-muscle-28432191.tsv",
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )


