perform_chi2_tests<- function(datasets, total_genes, verbose = FALSE) {
  n <- length(datasets)
  dataset_lengths <- sapply(datasets, length)
  names_list <- names(datasets)
  
  # Initialize matrices with names
  overlap_matrix <- matrix(0, n, n, dimnames = list(names_list, names_list))
  p_value_matrix <- matrix(0, n, n, dimnames = list(names_list, names_list))
  chi2_value_matrix <- matrix(0, n, n, dimnames = list(names_list, names_list))
  overlap_genes_matrix <- matrix("", n, n, dimnames = list(names_list, names_list))
  
  combinations <- combn(names_list, 2, simplify = FALSE)
  
  results <- lapply(combinations, function(pair) {
    i <- match(pair[1], names_list)
    j <- match(pair[2], names_list)
    
    overlapping_genes <- intersect(datasets[[i]], datasets[[j]])
    overlap_count <- length(overlapping_genes)
    external_genes_count <- length(total_genes) - length(unique(c(datasets[[i]], datasets[[j]])))
    
    matrix_chi2 <- matrix(c(overlap_count,
                            dataset_lengths[i] - overlap_count,
                            dataset_lengths[j] - overlap_count,
                            external_genes_count), nrow = 2)
    
    test_result <- chisq.test(matrix_chi2)
    
    if (verbose) {
      cat("Performed chi-square test on", pair[1], "and", pair[2], "\n")
    }
    
    list(
      i = i,
      j = j,
      overlap_count = overlap_count,
      p_value = test_result$p.value,
      chi2_value = test_result$statistic,
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
    chi2_value_matrix[i, j] <- result$chi2_value
    chi2_value_matrix[j, i] <- result$chi2_value
    overlap_genes_matrix[i, j] <- paste(result$overlapping_genes, collapse = ",")
    overlap_genes_matrix[j, i] <- overlap_genes_matrix[i, j]
  }
  
  return(list(
    number_overlap_matrix = overlap_matrix,
    p_value_matrix = p_value_matrix,
    chi2_value_matrix = chi2_value_matrix,
    overlap_genes_matrix = overlap_genes_matrix
  ))
}


# perform_chi2_tests_revised <- function(datasets, total_genes, verbose = FALSE) {
#   n <- length(datasets)
#   dataset_lengths <- sapply(datasets, length)
#   names_list <- names(datasets)
#   
#   # Initialize matrices with names
#   overlap_matrix <- matrix(0, n, n, dimnames = list(names_list, names_list))
#   p_value_matrix <- matrix(0, n, n, dimnames = list(names_list, names_list))
#   chi2_value_matrix <- matrix(0, n, n, dimnames = list(names_list, names_list))
#   overlap_genes_matrix <- matrix("", n, n, dimnames = list(names_list, names_list))
#   
#   for (i in 1:(n-1)) {
#     for (j in (i+1):n) {
#       overlapping_genes <- intersect(datasets[[i]], datasets[[j]])
#       overlap_count <- length(overlapping_genes)
#       overlap_matrix[i, j] <- overlap_count
#       overlap_matrix[j, i] <- overlap_count
#       overlap_genes_matrix[i, j] <- paste(overlapping_genes, collapse = ",")
#       overlap_genes_matrix[j, i] <- overlap_genes_matrix[i, j]
#       
#       external_genes_count <- length(total_genes) - length(unique(c(datasets[[i]], datasets[[j]])))
#       
#       matrix_chi2 <- matrix(c(overlap_count,
#                               dataset_lengths[i] - overlap_count,
#                               dataset_lengths[j] - overlap_count,
#                               external_genes_count), nrow = 2)
#       
#       test_result <- chisq.test(matrix_chi2)
#       p_value_matrix[i, j] <- test_result$p.value
#       p_value_matrix[j, i] <- test_result$p.value
#       chi2_value_matrix[i, j] <- test_result$statistic
#       chi2_value_matrix[j, i] <- test_result$statistic
#       
#       if (verbose) {
#         cat("Performed chi-square test on", names_list[i], "and", names_list[j], "\n")
#       }
#     }
#   }
#   
#   return(list(
#     number_overlap_matrix = overlap_matrix,
#     p_value_matrix = p_value_matrix,
#     chi2_value_matrix = chi2_value_matrix,
#     overlap_genes_matrix = overlap_genes_matrix
#   ))
# }
# 
# perform_chi2_tests <- function(datasets, total_genes, verbose = FALSE) {
# 
#   # Prepare matrices for overlap counts, p-values, chi2 values, and overlapping gene names
#   n <- length(datasets)
#   overlap_matrix <- matrix(0, n, n, dimnames = list(names(datasets), names(datasets)))
#   p_value_matrix <- matrix(0, n, n, dimnames = list(names(datasets), names(datasets)))
#   chi2_value_matrix <- matrix(0, n, n, dimnames = list(names(datasets), names(datasets)))
#   overlap_genes_matrix <- matrix("", n, n, dimnames = list(names(datasets), names(datasets)))
# 
#   for (i in 1:n) {
#     for (j in 1:n) {
#       overlapping_genes <- intersect(datasets[[i]], datasets[[j]])
#       overlap_matrix[i,j] <- length(overlapping_genes)
#       overlap_genes_matrix[i,j] <- paste(overlapping_genes, collapse = ",")
# 
#       intersect_genes <- unique(c(datasets[[i]], datasets[[j]]))
#       external_genes <- total_genes[!total_genes %in% intersect_genes]
# 
#       matrix_chi2 <- matrix(c(overlap_matrix[i,j],
#                               length(datasets[[i]]) - overlap_matrix[i,j],
#                               length(datasets[[j]]) - overlap_matrix[i,j],
#                               length(external_genes)), nrow = 2)
# 
#       test_result <- chisq.test(matrix_chi2)
#       p_value_matrix[i,j] <- test_result$p.value
#       chi2_value_matrix[i,j] <- test_result$statistic
# 
#       # Print out a message only if verbose is TRUE
#       if (verbose) {
#         cat("Performed chi-square test on", names(datasets)[i], "and", names(datasets)[j], "\n")
#       }
#     }
#   }
# 
#   return(list(
#     number_overlap_matrix = overlap_matrix,
#     p_value_matrix = p_value_matrix,
#     chi2_value_matrix = chi2_value_matrix,
#     overlap_genes_matrix = overlap_genes_matrix
#   ))
# }
