analyze_random_gene_sets <- function(reference_hgnc_vector, size_reference_list, comparison_gene_list, overlap_threshold, fdr_threshold, seed = NULL) {
  # Function: analyze_random_gene_sets
  # Description: 
  #   Performs an analysis on gene sets by generating random gene sets based on HGNC symbols,
  #   combining these with an additional gene set for comparison, and then performing chi-square tests 
  #   and overlap analysis. Useful for comparing randomly generated gene sets with other sets for 
  #   statistical and overlap analysis.
  # Parameters:
  #   reference_hgnc_vector - Vector of HGNC gene symbols.
  #   size_reference_list - List used to determine the sizes of the random gene sets.
  #   comparison_gene_list - Gene list used for additional comparison in the analysis.
  #   overlap_threshold - Threshold for determining significant overlap.
  #   fdr_threshold - False discovery rate threshold for statistical significance.
  #   seed - Optional; integer for setting the random number generator's seed for reproducibility.
  # Returns:
  #   A data structure containing the results of the chi-square tests and overlap analysis.
  
  # Generate random gene sets with the specified seed
  random_genes_list <- generate_random_gene_sets(
    genes_vector = reference_hgnc_vector,
    size_vector = sapply(size_reference_list, length),
    seed = seed
  )
  
  # Combine with comparison gene list and get unique values
  combined_genes_list <- c(random_genes_list, comparison_gene_list) %>% 
    lapply(., unique)
  
  # Perform chi-square tests
  chi2_results <- perform_chi2_tests(combined_genes_list, reference_hgnc_vector)
  
  # Process overlap results
  analysis_results <- processing_overlap_results(
    data = chi2_results,
    rows_to_filter = !rownames(chi2_results$p_value_matrix) %in% names(random_genes_list),
    cols_to_filter = names(random_genes_list),
    overlap_threshold = overlap_threshold,
    fdr_threshold = fdr_threshold,
    genes_list = combined_genes_list
  )
  
  rm(chi2_results, random_genes_list)
  
  gc()
  
  return(analysis_results)
}

