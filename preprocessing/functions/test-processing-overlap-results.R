processing_overlap_results <- function(data, rows_to_filter, cols_to_filter,
                                       fdr_threshold = 0.05, overlap_threshold = 2){
  # processing_overlap_results
  #
  # This function processes overlap results by filtering matrices based on provided row and column filters,
  # adjusting p-values, and organizing the results into a structured list containing original data,
  # significant data, and uniquely significant data.
  #
  # Parameters:
  #   data: A list of matrices where each matrix contains p-value data to be filtered and adjusted.
  #   rows_to_filter: A vector indicating which rows should be included after filtering.
  #   cols_to_filter: A vector indicating which columns should be included after filtering.
  #
  # Returns:
  #   A list containing three nested lists: 'original_data', 'significant_data', and 'significant_uniq_data'.
  #   Each nested list contains the filtered list of matrices, a data frame of the filtered results,
  #   and vectors of the filtered row and column names.

    
  # Define the function to filter matrices
  # This internal function takes a list of matrices and filters each matrix based on the specified rows and columns.
  filter_matrices <- function(data_list, rows_to_filter, cols_to_filter) {
    # Initialize an empty list to store the filtered matrices
    filtered_list_significant <- list()
    
    # Loop through each item in the matrix list by their names
    for (name in names(data_list)) {
      # Extract the matrix using its name
      matrix <- data_list[[name]]
      # Apply the filter using the provided indices for rows and columns
      filtered_matrix <- matrix[rows_to_filter, cols_to_filter]
      # Add the filtered matrix to the list with the same name
      filtered_list_significant[[name]] <- filtered_matrix
    }
    
    # Return the list of filtered matrices
    return(filtered_list_significant)
  }
  
  # Define the function to adjust p-values
  # This internal function takes a data list and applies False Discovery Rate (FDR) adjustment to the p-values.
  adjust_p_values <- function(data) {
    # Extract the p_value_matrix from the data
    p_value_matrix <- data$p_value_matrix
    
    # Apply the filters and adjustments
    fdr_value_matrix <- p_value_matrix %>%
      melt() %>%
      mutate(value = p.adjust(value, method = "fdr")) %>%
      dcast(Var1 ~ Var2, value.var = "value") %>%
      column_to_rownames(var = "Var1") %>%
      as.matrix()
    
    # Add the adjusted matrix back to the data list
    data$fdr_value_matrix <- fdr_value_matrix
    
    # Return the data list with the adjusted matrix
    return(data)
  }
  
  # The body of the main function starts here
  
  # 1. Prepare original_data by filtering and adjusting p-values
  # This section extracts data, applies row and column filters, and performs FDR adjustment.
  extract_data(
    data,
    rows_to_filter = rows_to_filter,
    cols_to_filter = cols_to_filter
  ) %>% 
    mutate(Var1 = as.character(Var1),
           Var2 = as.character(Var2)) %>% 
    mutate(fdr = p.adjust(p_value, method = "fdr")) -> original_df
  
  # Extract overlap genes
  original_df$overlap_genes %>% 
    strsplit(., split = ",") %>% 
    unlist %>% 
    unique() -> original_overlap_genes
  
  # Extract unique row and column names from the original data frame
  original_df$Var1 %>% unique() -> original_rows
  original_df$Var2 %>% unique() -> original_cols
  
  # Filter the original list of matrices and adjust p-values
  filter_matrices(data = data,
                  rows_to_filter = original_rows,
                  cols_to_filter = original_cols) %>% 
    adjust_p_values() -> original_list
  
  # 2. Prepare significant_data by filtering for significant results
  # This section filters the original data frame for significant results based on FDR and overlap criteria.
  original_df %>% 
    filter(fdr < fdr_threshold) %>% 
    filter(number_overlap >= overlap_threshold) -> significant_df
  
  # Extract overlap genes
  significant_df$overlap_genes %>% 
    strsplit(., split = ",") %>% 
    unlist %>% 
    unique() -> significant_overlap_genes
  
  # Extract unique row and column names from the significant data frame
  significant_df$Var1 %>% unique() -> significant_rows
  significant_df$Var2 %>% unique() -> significant_cols
  
  # Filter the list of matrices for significant results
  filter_matrices(data = original_list,
                  rows_to_filter = significant_rows,
                  cols_to_filter = original_cols) -> significant_list
  
  # 3. Prepare significant_uniq_data by selecting unique significant results
  # This section processes the significant data frame to select unique significant results.
  significant_df %>% 
    group_by(overlap_genes, Var2) %>% 
    nest() %>% 
    dplyr::mutate(data = map(data, ~ .x %>% 
                               arrange(fdr) %>% 
                               dplyr::slice(1))) %>% 
    unnest() %>% 
    ungroup %>% 
    as.data.frame() %>% 
    select(c(Var1, Var2, p_value, chi2, number_overlap, overlap_genes, fdr))-> significant_uniq_df
  
  # Extract overlap genes
  significant_uniq_df$overlap_genes %>% 
    strsplit(., split = ",") %>% 
    unlist %>% 
    unique() -> significant_uniq_overlap_genes
  
  # Extract unique row and column names from the uniquely significant data frame
  significant_uniq_df$Var1 %>% unique() -> significant_uniq_rows
  significant_uniq_df$Var2 %>% unique() -> significant_uniq_cols
  
  # Filter the original list of matrices for uniquely significant results
  filter_matrices(data = original_list,
                  rows_to_filter = significant_uniq_rows,
                  cols_to_filter = original_cols) -> significant_uniq_list
  
  # Prepare the output list containing all processed data
  output_list <- list(
    original_data = list(list = original_list, df = original_df, rows = original_rows, cols = original_rows, overlap_genes = original_overlap_genes),
    significant_data = list(list = significant_list, df = significant_df, rows = significant_rows, cols = significant_cols, overlap_genes = significant_overlap_genes),
    significant_uniq_data = list(list = significant_uniq_list, df = significant_uniq_df, rows = significant_uniq_rows, cols = significant_uniq_cols, overlap_genes = significant_uniq_overlap_genes)
  )
  
  # Return the structured list containing all processed data
  return(output_list)
}

processing_overlap_results(data = chi2_results_phenotypes,
                           rows_to_filter = !rownames(chi2_results_phenotypes$p_value_matrix) %in% tissues_clusters,
                           cols_to_filter = tissues_clusters[10:27]) -> tmp
# 
# tmp$original_data$list$fdr_value_matrix %>% .["biobankuk-warfarin-both_sexes--na-EUR-1e-08", "marpiech_tissues_dex_8"]

tmp$significant_uniq_data$overlap_genes



# Extract significant results for clusters vs papers
extract_data(
  chi2_results_phenotypes,
  !rownames(chi2_results_phenotypes$p_value_matrix) %in% tissues_clusters,
  tissues_clusters[10:27]
) %>% 
  dplyr::rename(phenotypes = Var1, cluster = Var2) %>% 
  mutate(fdr = p.adjust(p_value, method = "fdr")) %>%
  arrange(fdr) %>%
  filter(fdr < 0.05) %>% 
  filter(number_overlap > 1) %>% 
  group_by(overlap_genes, cluster) %>% 
  nest() %>% 
  dplyr::mutate(data = map(data, ~ .x %>% 
                             arrange(fdr) %>% 
                             dplyr::slice(1))) %>% 
  unnest() %>% 
  ungroup %>% 
  as.data.frame() %>% 
  .$phenotypes %>% 
  unique %>% 
  as.character() 