extract_keys_values <- function(data, column, keys) {
  # This function extracts specific key-value pairs from a specified column in a data frame.
  # The function splits the specified column by '|' to separate key-value pairs. 
  # Then, for each row, it finds the value corresponding to the provided key.
  # If the key doesn't exist in a row, NA is returned.
  # The function then adds these values as a new column in the data frame and returns the updated data frame.
  
  # Split the column by '|' to separate key-value pairs.
  split_info <- strsplit(as.character(data[[column]]), "\\|")
  
  # For each key, extract the corresponding values.
  for (key in keys) {
    
    # For each row, find the value corresponding to the provided key.
    values <- sapply(split_info, function(x) {
      
      # Split the key-value pairs by ':' to separate the keys and the values.
      key_value_pairs <- lapply(x, function(y) {
        if (grepl(":", y)) {
          strsplit(y, ":")[[1]]
        } else {
          return(c(y, NA))
        }
      })
      
      # Extract the keys and values.
      keys <- sapply(key_value_pairs, function(z) if(length(z) > 0) z[1] else NA)
      values <- sapply(key_value_pairs, function(z) if(length(z) > 1) z[2] else NA)
      
      # If the provided key exists, return the corresponding value.
      # If the key doesn't exist, return NA.
      if (key %in% keys) {
        return(values[keys == key])
      } else {
        return(NA)
      }
    })
    
    # Add the values as a new column in the data frame.
    data[[key]] <- unlist(values)
  }
  
  # Return the updated data frame.
  return(data)
}


print_simple_database <- function(data,
                                  remove_cols = c(
                                    'Info',
                                    'Ensembl_transcript_id',
                                    'Alias',
                                    'HGNC_symbol',
                                    'Ensembl_transcript_id',
                                    'RefSeq_mRNA_id'
                                  )) {
  
  # Function Description:
  # This function prints a simplified version of the database, 
  # removing specified columns for clarity.
  # By default, the function removes 'Info', 'Ensembl_transcript_id', and 'Alias' columns.
  
  # Check if all the columns specified exist in the data frame.
  if (!all(remove_cols %in% colnames(data))) {
    stop("Some columns specified for removal do not exist in the data frame.")
  }
  
  # Remove the specified columns from the data frame.
  data <- data[ , !(colnames(data) %in% remove_cols)]
  
  return(data)
}


pairwise_gene_overlap <- function(df, gene_col, label_col) {
  # Split the data frame into a list of gene vectors for each label
  gene_lists <- split(df[[gene_col]], df[[label_col]])
  
  # Get the names of all labels
  labels <- names(gene_lists)
  
  # Initialize a matrix to store the overlap counts
  overlap_matrix <- matrix(nrow = length(labels), ncol = length(labels))
  rownames(overlap_matrix) <- labels
  colnames(overlap_matrix) <- labels
  
  # Calculate the overlap for each pair of labels
  for (i in seq_along(gene_lists)) {
    for (j in seq_along(gene_lists)) {
      overlap_matrix[i, j] <- length(intersect(gene_lists[[i]], gene_lists[[j]]))
    }
  }
  
  return(overlap_matrix)
}


perform_chi2_tests <- function(datasets, total_genes) {
  # Prepare a square matrix to store the number of overlapping genes between datasets.
  # The dimension of the matrix is the square of the number of datasets. The row and 
  # column names of the matrix are the names of the datasets.
  overlap_matrix <- matrix(nrow = length(datasets), ncol = length(datasets),
                           dimnames = list(names(datasets), names(datasets)))
  
  # Calculate the number of overlapping genes for each pair of datasets. 
  # 'intersect()' function is used to find common elements between two datasets. 
  # The number of common elements is stored in the 'overlap_matrix'.
  for (i in 1:length(datasets)) {
    for (j in 1:length(datasets)) {
      overlap_matrix[i,j] <- length(intersect(datasets[[i]], datasets[[j]]))
    }
  }
  
  # Prepare a square matrix to store p-values of the chi-square test. 
  # The dimension and names are the same as the 'overlap_matrix'.
  p_value_matrix <- matrix(nrow = length(datasets), ncol = length(datasets),
                           dimnames = list(names(datasets), names(datasets)))
  
  # Prepare a square matrix to store the statistics of the chi-square test. 
  # The dimension and names are the same as the 'overlap_matrix'.
  chi2_value_matrix <- matrix(nrow = length(datasets), ncol = length(datasets),
                              dimnames = list(names(datasets), names(datasets)))
  
  # Perform chi-square test for each pair of datasets. 
  # For each pair, the number of intersecting genes, the number of genes unique to each dataset, 
  # and the number of genes outside both datasets are calculated. 
  # These numbers are used to create a 2x2 contingency table, which is the input of 'chisq.test()'.
  for (i in 1:length(datasets)) {
    for (j in 1:length(datasets)) {
      # Combine two datasets and extract unique genes
      intersect_genes <- unique(c(datasets[[i]], datasets[[j]]))  
      
      # Find genes in 'total_genes' that are not in 'intersect_genes'
      external_genes <- total_genes[!total_genes %in% intersect_genes]
      
      # Create a 2x2 contingency table. Each element of the table is the count of a specific category of genes.
      matrix_chi2 <- matrix(c(overlap_matrix[i,j], length(datasets[[i]]) - overlap_matrix[i,j], 
                              length(datasets[[j]]) - overlap_matrix[i,j],  length(external_genes)),
                            nrow = 2)
      
      # Perform chi-square test on the contingency table
      test_result <- chisq.test(matrix_chi2)
      
      # Store the p-value and the statistics of the chi-square test in their respective matrices
      p_value_matrix[i,j] <- test_result$p.value
      chi2_value_matrix[i,j] <- test_result$statistic 
      
      # Print out a message indicating the completion of a chi-square test
      cat("Performed chi-square test on", names(datasets)[i], "and", names(datasets)[j], "\n")
    }
  }
  
  # Return the matrices of overlapping gene numbers, p-values, and chi-square statistics
  return(list(number_overlap_matrix =  overlap_matrix, 
              p_value_matrix = p_value_matrix, 
              chi2_value_matrix = chi2_value_matrix))
}

perform_chi2_tests <- function(datasets, total_genes) {
  
  # Prepare matrices for overlap counts, p-values, chi2 values, and overlapping gene names
  n <- length(datasets)
  overlap_matrix <- matrix(0, n, n, dimnames = list(names(datasets), names(datasets)))
  p_value_matrix <- matrix(0, n, n, dimnames = list(names(datasets), names(datasets)))
  chi2_value_matrix <- matrix(0, n, n, dimnames = list(names(datasets), names(datasets)))
  overlap_genes_matrix <- matrix("", n, n, dimnames = list(names(datasets), names(datasets)))
  
  for (i in 1:n) {
    for (j in 1:n) {
      overlapping_genes <- intersect(datasets[[i]], datasets[[j]])
      overlap_matrix[i,j] <- length(overlapping_genes)
      overlap_genes_matrix[i,j] <- paste(overlapping_genes, collapse = ",")
      
      intersect_genes <- unique(c(datasets[[i]], datasets[[j]]))
      external_genes <- total_genes[!total_genes %in% intersect_genes]
      
      matrix_chi2 <- matrix(c(overlap_matrix[i,j], 
                              length(datasets[[i]]) - overlap_matrix[i,j], 
                              length(datasets[[j]]) - overlap_matrix[i,j],  
                              length(external_genes)), nrow = 2)
      
      test_result <- chisq.test(matrix_chi2)
      p_value_matrix[i,j] <- test_result$p.value
      chi2_value_matrix[i,j] <- test_result$statistic
      
      cat("Performed chi-square test on", names(datasets)[i], "and", names(datasets)[j], "\n")
    }
  }
  
  return(list(
    number_overlap_matrix = overlap_matrix, 
    p_value_matrix = p_value_matrix, 
    chi2_value_matrix = chi2_value_matrix,
    overlap_genes_matrix = overlap_genes_matrix
  ))
}


# Function to extract data from results_metabolism
extract_data <- function(results, rows_to_filter, cols_to_filter) {
  
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
  overlap_df <- melt_matrix(results$number_overlap_matrix, "number_overlap")
  
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

# function to filter chi2 results from list
# function to filter chi2 results from list
filter_matrices <- function(data, rows_to_filter, cols_to_filter, fdr_threshold = 0.05, overlap_threshold = 2) {
  
  data$p_value_matrix %>% 
    .[rows_to_filter, cols_to_filter] %>% 
    melt() %>% 
    mutate(value = p.adjust(value, method = "fdr")) %>% 
    dcast(Var1 ~ Var2, value.var = "value") %>%
    column_to_rownames(var = "Var1") %>% 
    as.matrix() -> fdr_value_matrix
  
  extract_data(
    data,
    rows_to_filter, cols_to_filter
  ) %>%
    mutate(fdr = p.adjust(p_value, method = "fdr")) %>%
    filter(fdr < fdr_threshold) %>% 
    filter(number_overlap >= overlap_threshold) %>%
    group_by(overlap_genes, Var2) %>% 
    nest() %>% 
    mutate(data = map(data, ~ .x %>% 
                        arrange(fdr) %>% 
                        slice(1))) %>% 
    unnest() %>% 
    ungroup %>% 
    as.data.frame() %>% 
    .$Var1 %>% 
    as.character %>% 
    unique() -> significant_res
  
  print(significant_res)
  
  
  # Initialize an empty list to store the filtered matrices
  filtered_list <- list()
  
  # Loop through each item in the matrix list by their names
  for (name in names(data)) {
    # Extract the matrix using its name
    matrix <- data[[name]]
    # Apply the filter using the provided indices for phenotypes and tissue clusters
    filtered_matrix <- matrix[rows_to_filter, cols_to_filter]
    # Add the filtered matrix to the list with the same name
    filtered_list[[name]] <- filtered_matrix
  }
  
  filtered_list$fdr_value_matrix <- fdr_value_matrix
  
  # Initialize an empty list to store the filtered matrices
  filtered_list_significant <- list()
  
  # Loop through each item in the matrix list by their names
  for (name in names(filtered_list)) {
    # Extract the matrix using its name
    matrix <- filtered_list[[name]]
    # Apply the filter using the provided indices for phenotypes and tissue clusters
    filtered_matrix <- matrix[significant_res, cols_to_filter]
    # Add the filtered matrix to the list with the same name
    filtered_list_significant[[name]] <- filtered_matrix
  }
  
  # Return the filtered list with the same names as the original list elements
  return(filtered_list_significant)
}


# This function takes a list of matrices and optionally two mapping vectors for column and row names.
replace_names_in_data <- function(data, col_mapping_vector=NULL, row_mapping_vector=NULL) {
  # This function takes a list of matrices and optionally two mapping vectors for column and row names.
  # It applies the mappings to the matrices to update their column and/or row names.
  # If a mapping vector is not provided or is NULL, the corresponding names (column or row) are not changed.
  # Parameters:
  #   data: A list of matrices whose names you want to update.
  #   col_mapping_vector: An optional named vector where names are old column names and values are new column names.
  #   row_mapping_vector: An optional named vector where names are old row names and values are new row names.
  # Returns:
  #   A list of matrices with updated column and/or row names based on the provided mappings.
  
  
  # Inner function to replace column and row names for a single matrix
  replace_names <- function(matrix, col_mapping_vector, row_mapping_vector) {
    # Replace column names if a column mapping vector is provided
    if (!is.null(col_mapping_vector) && length(col_mapping_vector) > 0) {
      cols_to_replace <- colnames(matrix) %in% names(col_mapping_vector)
      new_colnames <- colnames(matrix)
      new_colnames[cols_to_replace] <- col_mapping_vector[colnames(matrix)[cols_to_replace]]
      colnames(matrix) <- new_colnames
    }
    
    # Replace row names if a row mapping vector is provided
    if (!is.null(row_mapping_vector) && length(row_mapping_vector) > 0) {
      rows_to_replace <- rownames(matrix) %in% names(row_mapping_vector)
      new_rownames <- rownames(matrix)
      new_rownames[rows_to_replace] <- row_mapping_vector[rownames(matrix)[rows_to_replace]]
      rownames(matrix) <- new_rownames
    }
    
    # Return the matrix with updated names
    return(matrix)
  }
  
  # Apply the replace_names function to each matrix in the data list
  updated_data <- lapply(data, replace_names, col_mapping_vector=col_mapping_vector, row_mapping_vector=row_mapping_vector)
  
  # Return the list of matrices with updated names
  return(updated_data)
}

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
