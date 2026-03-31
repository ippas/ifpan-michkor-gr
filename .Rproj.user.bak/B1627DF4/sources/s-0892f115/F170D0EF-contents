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


