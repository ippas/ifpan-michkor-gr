perform_overlap_permutation_analysis_multicore <- function(permutations,
                                                           seed = NULL,
                                                           num_cores = 1,
                                                           reference_hgnc_vector,
                                                           size_reference_list,
                                                           comparison_gene_list,
                                                           overlap_threshold,
                                                           fdr_threshold) {
  # Function: perform_overlap_permutation_analysis_multicore
  # Description: 
  #   Performs a parallel permutation analysis on gene sets using multicore processing.
  #   This function is designed for statistical and overlap analysis in genetic data,
  #   executing multiple iterations of gene set analysis across multiple cores.
  
  # Set up seeds for each permutation if a seed is provided
  if (!is.null(seed)) {
    set.seed(seed)
    permutation_seeds <- sample(1000000, size = permutations)
  } else {
    permutation_seeds <- rep(NULL, permutations)
  }
  
  # Define the function for a single permutation
  perform_single_permutation <- function(i, seed, reference_hgnc_vector, size_reference_list, 
                                         comparison_gene_list, overlap_threshold, fdr_threshold) {

    # Use tryCatch to handle errors
    tryCatch({
      # Perform random gene set analysis
      results <- analyze_random_gene_sets(
        seed = seed[i],
        reference_hgnc_vector = reference_hgnc_vector,
        size_reference_list = size_reference_list,
        comparison_gene_list = comparison_gene_list,
        overlap_threshold = overlap_threshold,
        fdr_threshold = fdr_threshold
      )

      # Return the number of rows
      nrow(results$significant_uniq_data$df)
    }, error = function(e) {
      # In case of an error, return 0
      return(0)
    })
  }
  
  # Execute permutations in parallel using mclapply
  permutation_results <- parallel::mclapply(1:permutations, perform_single_permutation, 
                                            seed = permutation_seeds, 
                                            reference_hgnc_vector = reference_hgnc_vector, 
                                            size_reference_list = size_reference_list, 
                                            comparison_gene_list = comparison_gene_list, 
                                            overlap_threshold = overlap_threshold, 
                                            fdr_threshold = fdr_threshold, 
                                            mc.cores = num_cores)
  
  return(permutation_results)
}

