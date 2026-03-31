perform_fisher_tests <- function(datasets, total_genes, verbose = FALSE) {
  n <- length(datasets)
  dataset_lengths <- sapply(datasets, length)
  names_list <- names(datasets)
  
  # Initialize matrices with names
  overlap_matrix <- matrix(0, n, n, dimnames = list(names_list, names_list))
  p_value_matrix <- matrix(0, n, n, dimnames = list(names_list, names_list))
  fisher_value_matrix <- matrix(0, n, n, dimnames = list(names_list, names_list))
  overlap_genes_matrix <- matrix("", n, n, dimnames = list(names_list, names_list))
  
  combinations <- combn(names_list, 2, simplify = FALSE)
  
  results <- lapply(combinations, function(pair) {
    i <- match(pair[1], names_list)
    j <- match(pair[2], names_list)
    
    overlapping_genes <- intersect(datasets[[i]], datasets[[j]])
    overlap_count <- length(overlapping_genes)
    external_genes_count <- length(total_genes) - length(unique(c(datasets[[i]], datasets[[j]])))
    
    matrix_test <- matrix(c(overlap_count,
                            dataset_lengths[i] - overlap_count,
                            dataset_lengths[j] - overlap_count,
                            external_genes_count), nrow = 2)
    
    test_result <- fisher.test(matrix_test)
    
    if (verbose) {
      cat("Performed Fisher's exact test on", pair[1], "and", pair[2], "\n")
    }
    
    list(
      i = i,
      j = j,
      overlap_count = overlap_count,
      p_value = test_result$p.value,
      fisher_value = test_result$estimate,
      overlapping_genes = overlapping_genes
    )
  })
  
  for (result in results) {
    i <- result$i
    j <- result$j
    overlap_matrix[i, j] <- result$overlap_count
    overlap_matrix[j, i] <- result$overlap_count
    p_value_matrix[i, j] <- result$p_value
    p_value_matrix[j, i] <- result$p_value
    fisher_value_matrix[i, j] <- result$fisher_value
    fisher_value_matrix[j, i] <- result$fisher_value
    overlap_genes_matrix[i, j] <- paste(result$overlapping_genes, collapse = ",")
    overlap_genes_matrix[j, i] <- overlap_genes_matrix[i, j]
  }
  
  return(list(
    number_overlap_matrix = overlap_matrix,
    p_value_matrix = p_value_matrix,
    fisher_value_matrix = fisher_value_matrix,
    overlap_genes_matrix = overlap_genes_matrix
  ))
}
