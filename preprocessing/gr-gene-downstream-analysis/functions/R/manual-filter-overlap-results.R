# firsr version of function
manual_filter_overlap_results <- function(data, rows_to_remove = c("pmid:24777604_embryos_hypothalamic-region_NPSCs_up",
                                                                  "pmid:24777604_embryos_hypothalamic-region_NPSCs-KO-Cav1_up",
                                                                  "pmid:26606517_embryos_hypothalamic-region_NPSCs_up" )) {
  # All row names
  all_rows <- data$significant_uniq_data$rows
  
  # Rows to keep
  rows_to_keep <- all_rows[!all_rows %in% rows_to_remove]
  
  data$significant_uniq_data$list$number_overlap_matrix[rows_to_keep,] -> number_overlap_matrix
  data$significant_uniq_data$list$p_value_matrix[rows_to_keep,] -> p_value_matrix
  data$significant_uniq_data$list$chi2_value_matrix[rows_to_keep,] -> chi2_value_matrix
  data$significant_uniq_data$list$overlap_genes_matrix[rows_to_keep,] -> overlap_genes_matrix
  data$significant_uniq_data$list$fdr_value_matrix[rows_to_keep,] -> fdr_value_matrix
  
  list <- list(
    number_overlap_matrix = number_overlap_matrix,
    p_value_matrix = p_value_matrix,
    chi2_value_matrix = chi2_value_matrix,
    overlap_genes_matrix = overlap_genes_matrix,
    fdr_value_matrix = fdr_value_matrix
  )
  
  df <- data$significant_uniq_data$df %>% filter(Var1 %in% rows_to_keep)
  rows <- data$significant_uniq_data$rows
  cols <- data$significant_uniq_data$cols
  overlap_genes <- df$overlap_genes %>% strsplit(., ",") %>% unlist %>% unique

  manual_filter <- list(df = df,
                        list = list,
                        rows = rows, 
                        cols = cols,
                        overlap_genes)
  
  data$manual_filter_overlap_results <- manual_filter 
  
  return(data)
} 

# second version
manual_filter_overlap_results <- function(data, rows_to_remove) {
  # Extracting the 'significant_uniq_data' component from the input data
  sig_uniq_data <- data$significant_uniq_data
  
  # Identifying rows to keep by excluding the rows specified in 'rows_to_remove'
  rows_to_keep <- setdiff(sig_uniq_data$rows, rows_to_remove)
  
  # Filtering matrices and dataframe based on 'rows_to_keep'
  filtered_list <- lapply(sig_uniq_data$list, function(matrix) matrix[rows_to_keep, ])
  filtered_df <- dplyr::filter(sig_uniq_data$df, Var1 %in% rows_to_keep)
  
  # Extracting unique overlap genes from the filtered dataframe
  overlap_genes <- unique(unlist(strsplit(filtered_df$overlap_genes, ",")))
  
  # Creating the manual filter list with updated components
  manual_filter <- list(
    df = filtered_df,
    list = filtered_list,
    rows = rows_to_keep, 
    cols = sig_uniq_data$cols,
    overlap_genes = overlap_genes
  )
  
  # Adding the manual filter list to the original data
  data$manual_filter_overlap_results <- manual_filter 
  
  # Returning the updated data
  return(data)
} 
 
# manual_filter_overlap_results(data = clusters_papers_data) -> tmp
# tmp$manual_filter_overlap_results$list$number_overlap_matrix
#           
# 
# clusters_papers_data$significant_uniq_data$df$overlap_genes %>% strsplit(., ",") %>% unlist %>% unique
