generate_random_gene_sets_with_frequency <- function(gene_list, genes_vector, seed = NULL) {
  # Set seed for reproducibility if specified
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Flatten the list and calculate frequencies
  all_genes <- unlist(gene_list)
  
  # Create a named vector of random genes directly, avoiding the need to calculate frequencies
  unique_genes <- unique(all_genes)
  random_genes <- setNames(sample(genes_vector, length(unique_genes), replace = TRUE), unique_genes)
  
  # Use lapply to iterate over the list and replace each gene with its random counterpart
  result_lists <- lapply(gene_list, function(gene_list_inner) {
    sapply(gene_list_inner, function(gene) random_genes[gene])
  })
  
  # Naming the result lists (optional, for consistency with your original function)
  names(result_lists) <- paste0("random", seq_along(result_lists))
  
  result_lists <- lapply(result_lists, unname)
  
  return(result_lists)
}
