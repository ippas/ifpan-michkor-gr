perform_overlap_permutation_analysis_parallel <- function(permutations,
                                                        seed = NULL,
                                                        num_cores = 1,
                                                        reference_hgnc_vector,
                                                        size_reference_list,
                                                        comparison_gene_list,
                                                        overlap_threshold,
                                                        fdr_threshold) {
  # Function: perform_overlap_permutation_analysis_future
  # Description: 
  #   Performs a parallel permutation analysis on gene sets using the 'future' package.
  #   This function is designed for statistical and overlap analysis in genetic data,
  #   executing multiple iterations of gene set analysis across multiple cores.
  # Parameters:
  #   permutations - Number of permutations to perform.
  #   seed - Optional seed for random number generation for reproducibility.
  #   num_cores - Number of cores to use for parallel processing.
  #   reference_hgnc_vector - Vector of HGNC gene symbols used as a reference.
  #   size_reference_list - List to determine sizes of random gene sets.
  #   comparison_gene_list - Gene list used for additional comparison.
  #   overlap_threshold - Threshold for determining significant overlap.
  #   fdr_threshold - False discovery rate threshold for statistical significance.
  # Returns:
  #   A list of results from each permutation, each representing the outcome of the gene set analysis.
  
  # Initialize a future plan for parallel processing
  # Initialize a future plan for parallel processing
  plan(future::multisession, workers = num_cores)
  # plan(future::multicore, workers = num_cores)
  # plan(future::multiprocess, workers = num_cores)
  # plan(workers = num_cores)
  
  # Create a function that performs a single permutation
  perform_single_permutation <- function(i) {
    # Existing code for setting seed and generating permutation seeds
    if (!is.null(seed)) {
      set.seed(seed)
      permutation_seeds <- sample(1000000, size = permutations)
    } else {
      permutation_seeds <- rep(NULL, permutations)
    }
    
    # Use tryCatch to handle errors
    tryCatch({
      # Perform random gene set analysis
      analyze_random_gene_sets(
        seed = permutation_seeds[i],
        reference_hgnc_vector = reference_hgnc_vector,
        size_reference_list = size_reference_list,
        comparison_gene_list = comparison_gene_list,
        overlap_threshold = overlap_threshold,
        fdr_threshold = fdr_threshold
      ) -> results
      
      # print(i)
      
      # Return the number of rows
      number_significant_results <- nrow(results$significant_uniq_data$df)
      rm(results)
      return(number_significant_results)
    }, error = function(e) {
      # In case of an error, return or print 0
      # Uncomment the line below to print the error message, if needed
      # print(e)
      return(0)
    })
  }
  
  # Execute permutations in parallel using future.apply's future_lapply
  permutation_results <- future.apply::future_lapply(1:permutations, perform_single_permutation)
  
  # Stop the future plan
  future::plan(NULL)
  
  return(permutation_results)
}
