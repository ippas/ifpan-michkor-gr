# Define the function with an additional `tissue_column` argument
categorize_tissue <- function(df, tissue_column) {
  # Check if the specified tissue_column exists in the dataframe
  if(!tissue_column %in% names(df)) {
    stop(paste("The dataframe does not contain a", tissue_column, "column."), call. = FALSE)
  }
  
  # Perform the mutations using the specified column for tissue information
  df <- df %>%
    mutate(system = case_when(
      str_detect(!!sym(tissue_column), "brain|cortex(?<!adrenal-cortex)|hippocampal|neuronal|astrocyte|striatum|hypothalamus|hippocampus|ARC|pituitary") ~ "brain_and_nervous_system",
      str_detect(!!sym(tissue_column), "blood|placenta") ~ "blood_and_immune_system",
      str_detect(!!sym(tissue_column), "lung") ~ "respiratory_system",
      str_detect(!!sym(tissue_column), "liver|kidney|small-intestine|spleen") ~ "digestive_system",
      str_detect(!!sym(tissue_column), "muscle|adipose|adipocyte|bone|cartilage|anterior") ~ "musculoskeletal_system",
      str_detect(!!sym(tissue_column), "embryos") ~ "embryos",
      str_detect(!!sym(tissue_column), "adrenal-gland|adrenal-cortex") ~ "endocrine_system",
      TRUE ~ "other_tissues"
    )) %>%
    mutate(simple_tissue = case_when(
      str_detect(!!sym(tissue_column), "brain|cortex(?<!adrenal-cortex)|hippocampal|neuronal|astrocyte|striatum|hypothalamus|hippocampus|ARC|pituitary") ~ "brain",
      grepl("embryos", !!sym(tissue_column), ignore.case = TRUE) ~ "embryos",
      grepl("placenta", !!sym(tissue_column), ignore.case = TRUE) ~ "placenta",
      grepl("bone", !!sym(tissue_column), ignore.case = TRUE) ~ "bone",
      grepl("cartilage", !!sym(tissue_column), ignore.case = TRUE) ~ "cartilage",
      grepl("liver", !!sym(tissue_column), ignore.case = TRUE) ~ "liver",
      grepl("kidney", !!sym(tissue_column), ignore.case = TRUE) ~ "kidney",
      grepl("spleen", !!sym(tissue_column), ignore.case = TRUE) ~ "spleen",
      grepl("adipose|adipocyte", !!sym(tissue_column), ignore.case = TRUE) ~ "adipose",
      grepl("blood", !!sym(tissue_column), ignore.case = TRUE) ~ "blood",
      grepl("lung", !!sym(tissue_column), ignore.case = TRUE) ~ "lung",
      grepl("adrenal-gland|adrenal-cortex", !!sym(tissue_column), ignore.case = TRUE) ~ "adrenal-gland",
      grepl("small-intestine", !!sym(tissue_column), ignore.case = TRUE) ~ "small-intestine",
      grepl("muscle|anterior", !!sym(tissue_column), ignore.case = TRUE) ~ "muscle",
      TRUE ~ "Other"  # This will categorize any tissue not matched above as "Other"
    )) 
  
  return(df)
}


################################################################################
process_tissue_data <- function(data, columns, tissue_column = "tissue", additional_columns = NULL) {
  library(dplyr)
  
  data_processed <- data %>%
    mutate(!!tissue_column := ifelse(is.na(!!sym(tissue_column)), cell, !!sym(tissue_column))) %>%
    categorize_tissue(tissue_column = tissue_column) %>%
    select(c(columns, additional_columns)) %>%
    unique()
  
  return(data_processed)
}

barplot_tissue_data <- function(data, count_column = "simple_tissue", threshold = 3) {
  library(ggplot2)
  
  data %>%
    .[[count_column]] %>%
    table() %>%
    as.data.frame() %>% 
    set_colnames(c("tissue", "freq")) %>% 
    arrange(desc(freq)) %>%
    ggplot(aes(x = reorder(tissue, -freq), y = freq)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = threshold, color = "red", linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}


process_tissue_data(
  data = gr_database_blocked_genes_list$GR_dependent_gene_lists_pmid_tissue_cell_treatment_type_time_dose_regulation$metadata,
  columns = c("source", "tissue", "system", "simple_tissue"),
  additional_columns = c("time") # Add or remove as needed
) %>% #select(c(source, simple_tissue)) %>% unique() %>% 
  barplot_tissue_data(count_column = "source")
  
