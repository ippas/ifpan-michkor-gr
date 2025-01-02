analyze_gene_list_overlap(
  col_lists = categorized_gene_lists$metabolome_gene_lists$gene_lists,
  row_lists = c(gr_database_blocked_gene_lists$master_and_tissue_gene_lists_down, gr_database_blocked_gene_lists$master_and_tissue_gene_lists_up),
  reference_hgnc_vector = hgnc_symbols_vector_v110, 
  keep_original_data = TRUE, 
  fdr_threshold = 0.01, 
  overlap_threshold = 2
) -> results_overlap_secondary_gene_lists$metabolon_vs_secondary_gene_lists

results_overlap_secondary_gene_lists$metabolon_vs_secondary_gene_lists$significant_data$df

combined_list <- lapply(c(categorized_gene_lists$metabolome_gene_lists$gene_lists, c(gr_database_blocked_gene_lists$master_and_tissue_gene_lists_down, gr_database_blocked_gene_lists$master_and_tissue_gene_lists_up)), unique)
col_lists = categorized_gene_lists$metabolome_gene_lists$gene_lists
row_lists = c(gr_database_blocked_gene_lists$master_and_tissue_gene_lists_down, gr_database_blocked_gene_lists$master_and_tissue_gene_lists_up)

analysis_results <-
  perform_chi2_tests(combined_list, hgnc_symbols_vector_v110)
rows_to_filter = !rownames(analysis_results$p_value_matrix) %in% names(row_lists)
cols_to_filter = names(row_lists)

melt_matrix <- function(matrix, column_name) {
  matrix %>%
    # Filter the matrix based on the provided rows and columns
    .[rows_to_filter, cols_to_filter] %>%
    # Melt the matrix to transform it into a long format
    melt() %>%
    # Rename the columns of the melted dataframe
    `colnames<-`(c("Var1", "Var2", column_name))
}

results <- analysis_results

p_value_df <- melt_matrix(results$p_value_matrix, "p_value")

# Extract chi2 values and store in a dataframe
chi2_df <- melt_matrix(results$chi2_value_matrix, "chi2")



# Extract number of overlaps and store in a dataframe
overlap_df <- melt_matrix(results$number_overlap_matrix, "gene_overlap_count")

# Extract overlapping gene names and store in a dataframe
overlap_genes_df <- melt_matrix(results$overlap_genes_matrix, "overlap_genes")

overlap_df %>% filter(gene_overlap_count > 1)

p_value_df %>%
  left_join(chi2_df, by = c("Var1", "Var2")) %>%dim
  left_join(overlap_df, by = c("Var1", "Var2")) %>% dim


p_value_df %>% head
chi2_df %>%  head
overlap_df %>% head

p_value_df %>% tail
chi2_df %>%  tail
overlap_df %>% tail
overlap_genes_df %>% dim

final_df <-  bind_cols(p_value_df, select(chi2_df, chi2), select(overlap_df, gene_overlap_count), select(overlap_genes_df, overlap_genes))

# Merge p_value_df and chi2_df
final_df <- merge(p_value_df, chi2_df, by = c("Var1", "Var2"))

final_df %>% dim

# Merge the above result with overlap_df
final_df <- merge(final_df, overlap_df, by = c("Var1", "Var2"))

final_df %>% dim

# Merge the above result with overlap_genes_df
final_df <- merge(final_df, overlap_genes_df, by = c("Var1", "Var2"))

final_df %>% dim
# funckajca extract_data -> tutaj jest jakiś problem, wynik tej fukcji daje jakieś dizwne wyniki, przy łączeniu final_df robi dziwne rzeczy
extract_data
function(results, rows_to_filter, cols_to_filter) {
  
  # Helper function to melt a given matrix and set column names
  # Args:
  #   matrix: The matrix to be melted
  #   column_name: The name for the value column after melting
  # Returns:
  #   A melted dataframe with columns "Var1", "Var2", and the specified column_name
  melt_matrix <- function(matrix, column_name) {
    matrix %>%
      # Filter the matrix based on the provided rows and columns
      .[rows_to_filter, cols_to_filter] %>%
      # Melt the matrix to transform it into a long format
      melt() %>%
      # Rename the columns of the melted dataframe
      `colnames<-`(c("Var1", "Var2", column_name))
  }
  
  # Extract data from each matrix in the results list using the melt_matrix function
  
  # Extract p-values and store in a dataframe
  p_value_df <- melt_matrix(results$p_value_matrix, "p_value")
  
  # Extract chi2 values and store in a dataframe
  chi2_df <- melt_matrix(results$chi2_value_matrix, "chi2")
  
  # Extract number of overlaps and store in a dataframe
  overlap_df <- melt_matrix(results$number_overlap_matrix, "gene_overlap_count")
  
  # Extract overlapping gene names and store in a dataframe
  overlap_genes_df <- melt_matrix(results$overlap_genes_matrix, "overlap_genes")
  
  # Merge the extracted dataframes by the "Var1" and "Var2" columns
  # This ensures that the final dataframe has a row for each pair of datasets
  # and columns for p_value, chi2, number_overlap, and overlap_genes
  
  # Merge p_value_df and chi2_df
  final_df <- merge(p_value_df, chi2_df, by = c("Var1", "Var2"))
  
  # Merge the above result with overlap_df
  final_df <- merge(final_df, overlap_df, by = c("Var1", "Var2"))
  
  # Merge the above result with overlap_genes_df
  final_df <- merge(final_df, overlap_genes_df, by = c("Var1", "Var2"))
  
  # Return the final merged dataframe
  return(final_df)
}

extract_data(
  analysis_results,
  rows_to_filter =!rownames(analysis_results$p_value_matrix) %in% names(row_lists),
  cols_to_filter = names(row_lists)
) %>% 
  mutate(Var1 = as.character(Var1),
         Var2 = as.character(Var2)) %>% filter(gene_overlap_count > 1)
  mutate(fdr = p.adjust(p_value, method = "fdr")) %>%  
  mutate(
    overlap_genes = overlap_genes %>% 
      strsplit(., ',') %>% 
      map(~sort(.) %>% paste(collapse = ',')) %>% unlist
  ) -> original_df


function(data, genes_list, rows_to_filter, cols_to_filter,
         fdr_threshold = 0.05, overlap_threshold = 2){


  
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
    mutate(fdr = p.adjust(p_value, method = "fdr")) %>% 
    mutate(
      overlap_genes = overlap_genes %>% 
        strsplit(., ',') %>% 
        map(~sort(.) %>% paste(collapse = ',')) %>% unlist
    ) -> original_df
  
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
    filter(gene_overlap_count >= overlap_threshold) -> significant_df
  
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
    select(c(Var1, Var2, p_value, chi2, gene_overlap_count, overlap_genes, fdr))-> significant_uniq_df
  
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
  
  gene_list_sizes <- genes_list %>% sapply(., length)
  
  # Prepare the output list containing all processed data
  output_list <- list(
    original_data = list(list = original_list, df = original_df, rows = original_rows, cols = original_rows, overlap_genes = original_overlap_genes),
    significant_data = list(list = significant_list, df = significant_df, rows = significant_rows, cols = significant_cols, overlap_genes = significant_overlap_genes),
    significant_uniq_data = list(list = significant_uniq_list, df = significant_uniq_df, rows = significant_uniq_rows, cols = significant_uniq_cols, overlap_genes = significant_uniq_overlap_genes),
    gene_list_sizes = gene_list_sizes
  )
  
  # Return the structured list containing all processed data
  return(output_list)
}





function(row_lists,
         col_lists,
         reference_hgnc_vector,
         fdr_threshold = 0.05,
         overlap_threshold = 3,
         keep_original_data = TRUE) {
  # Combine and get unique elements from both lists
  combined_list <- lapply(c(row_lists, col_lists), unique)
  
  # Perform chi2 tests or similar analysis
  analysis_results <-
    perform_chi2_tests(combined_list, reference_hgnc_vector)
  
  # Process overlap results or similar post-analysis processing
  processed_data <- processing_overlap_results(
    data = analysis_results,
    rows_to_filter = !rownames(analysis_results$p_value_matrix) %in% names(row_lists),
    cols_to_filter = names(row_lists),
    overlap_threshold = {{overlap_threshold}},
    fdr_threshold = {{fdr_threshold}},
    genes_list = combined_list
  )
  
  
  # Remove the "original_data" element from the list
  if(!keep_original_data){
    processed_data <- processed_data[!(names(processed_data) == "original_data")]
  }
  
  return(processed_data)
}

