perform_chi2_tests <- function(datasets, total_genes, verbose = FALSE) {

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

      # Print out a message only if verbose is TRUE
      if (verbose) {
        cat("Performed chi-square test on", names(datasets)[i], "and", names(datasets)[j], "\n")
      }
    }
  }

  return(list(
    number_overlap_matrix = overlap_matrix,
    p_value_matrix = p_value_matrix,
    chi2_value_matrix = chi2_value_matrix,
    overlap_genes_matrix = overlap_genes_matrix
  ))
}
