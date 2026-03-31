overlap_coefficient <- function(set1, set2) {
  if (length(set1) == 0 & length(set2) == 0) {
    return(1)  # if both sets are empty, define their similarity as 1
  }
  intersection_size <- length(intersect(set1, set2))
  smaller_set_size <- min(length(set1), length(set2))
  return(intersection_size / smaller_set_size)
}


overlap_coefficient_matrix <- function(gene_list) {
  num_genes <- length(gene_list)
  matrix_result <- matrix(0, nrow = num_genes, ncol = num_genes)
  for (i in 1:num_genes) {
    for (j in i:num_genes) {
      if (i == j) {
        matrix_result[i, j] <-NA
      } else {
        index <- overlap_coefficient(gene_list[[i]], gene_list[[j]])
        matrix_result[i, j] <- index
        matrix_result[j, i] <- index  # since the matrix is symmetric
      }
    }
  }
  colnames(matrix_result) <- rownames(matrix_result) <- names(gene_list)
  return(matrix_result)
}
