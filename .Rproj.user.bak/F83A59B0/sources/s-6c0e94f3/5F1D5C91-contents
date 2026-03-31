generate_random_gene_list_dependent <- function(genes_list, genes_vector, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  original_genes <- unique(unlist(unname(genes_list)))
  random_genes <- sample(genes_vector, length(original_genes))
  
  # Create a named vector for faster lookup
  gene_replacement <- setNames(random_genes, original_genes)
  
  # Replace genes using vectorized operations where possible
  replace_genes <- function(genes_vector, replacement_vector) {
    # Use the replacement vector for fast lookup
    replaced_genes <- replacement_vector[genes_vector]
    # Replace NA values with original genes (for genes not found in the mapping)
    replaced_genes[is.na(replaced_genes)] <- genes_vector[is.na(replaced_genes)]
    return(replaced_genes)
  }
  
  random_gene_list <- lapply(genes_list, replace_genes, replacement_vector = gene_replacement)
  
  names(random_gene_list) <- paste0("random", seq_along(random_gene_list))
  
  return(random_gene_list)
}
