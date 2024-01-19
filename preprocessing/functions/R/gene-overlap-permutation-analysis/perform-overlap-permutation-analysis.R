perform_overlap_permutation_analysis <- function(permutations, seed = NULL, ...) {
  # Function: perform_overlap_permutation_analysis
  # Description: 
  #   This function performs a permutation analysis on gene sets. It is designed to execute
  #   a specified number of permutations, each involving a random gene set analysis. The function
  #   is particularly useful for statistical and overlap analysis in genetic data, where multiple
  #   iterations of the analysis are required to assess variability and significance.
  # Parameters:
  #   permutations - The number of permutations to perform. Each permutation involves a separate
  #                  execution of the gene set analysis.
  #   seed - Optional; an integer value used to set the seed for random number generation, ensuring
  #          reproducibility of results. If NULL, random number generation is not controlled by a fixed seed.
  #   ... - Additional arguments to be passed to the analyze_random_gene_sets function.
  # Returns:
  #   A list of results from each permutation. Each element in the list corresponds to the result
  #   of one permutation of the gene set analysis.
  
  # Initialize a list to store results of each permutation
  permutation_results <- list()
  
  # Set initial seed for reproducibility
  if (!is.null(seed)) {
    set.seed(seed)
    permutation_seeds <- sample(1000000, size = permutations)
  } else {
    permutation_seeds <- rep(NULL, permutations)
  }
  
  # Execute permutations
  for (i in 1:permutations) {
    # Generate a seed for each permutation if initial seed is provided
    # permutation_seed <- if (!is.null(seed)) permutation_seed[i] else NULL
    
    # Perform random gene set analysis with additional arguments passed via ...
    result <- analyze_random_gene_sets(
      seed = permutation_seed[i],
      ...
    )
    
    # Store the result
    permutation_results[[i]] <- result
  }
  
  return(permutation_results)
}