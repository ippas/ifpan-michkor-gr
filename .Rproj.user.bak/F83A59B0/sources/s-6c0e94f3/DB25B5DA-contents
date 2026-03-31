generate_random_gene_list_independent <- function(genes_vector, size_vector, seed = NULL) {
  # Function: generate_random_gene_sets
  # Description: 
  #   Generates a list of random gene sets. Each set is a vector of gene names randomly 
  #   selected from a provided vector of genes. Allows specification of the size of each 
  #   gene set and optionally takes a seed for random number generation for reproducibility.
  # Parameters:
  #   genes_vector - Vector containing gene names for sampling.
  #   size_vector - Vector of integers specifying the size of each random gene set.
  #   seed - Optional; integer for setting the seed for random number generation.
  # Returns:
  #   A named list of gene vectors, each representing a random gene set. 
  #   Names are formatted as "random1", "random2", etc.
  
  # Set seed for reproducible random generation if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Generate random gene sets
  random_gene_sets <- lapply(seq_along(size_vector), function(i) {
    set_size <- size_vector[i]
    genes <- sample(genes_vector, set_size)
    genes
  })
  
  # Assign names to each gene set
  names(random_gene_sets) <- paste("random", seq_along(random_gene_sets), sep = "")
  return(random_gene_sets)
}