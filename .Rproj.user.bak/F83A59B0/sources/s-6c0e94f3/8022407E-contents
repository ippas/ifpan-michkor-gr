analyze_uk_biobank_categories_gene_overlap <- function(categories, genes_phenotypes, papers_gene_list, hgnc_symbols_vector, path_metafile, fdr_threshold = 0.01, overlap_threshold = 3, permutations = 1000, seed = 123, num_cores = 30, random_genes_type = "dependent") {
  # Initialize an empty list to store results
  overlap_results <- list()
  
  # Start time measurement
  start_time <- Sys.time()
  
  # Sequential processing using a for loop
  for (category in categories) {
    
    print(category)
    loop_start <- Sys.time()
    
    # Filter genes_phenotypes based on the category
    col_list <- filter_phenotypes_by_category(
      genes_list = genes_phenotypes,
      path_metafile = path_metafile,
      category_name = category
    )
    
    # Prepare row_list by excluding specific entries from papers_gene_list
    row_list <- papers_gene_list[!grepl("marpiech", names(papers_gene_list))]
    
    # Use tryCatch to handle errors during gene list overlap analysis
    category_overlap_results <- tryCatch({
      analyze_gene_list_overlap(
        row_lists = row_list,
        col_lists = col_list,
        reference_hgnc_vector = hgnc_symbols_vector,
        fdr_threshold = fdr_threshold,
        overlap_threshold = overlap_threshold,
        keep_original_data = TRUE
      )
    }, error = function(e) {
      warning(sprintf("Error in analyze_gene_list_overlap for category '%s': %s", category, e$message))
      return(NULL)  # Return NULL in case of error
    })
    
    # Skip further processing if category_overlap_results is NULL
    if (is.null(category_overlap_results)) {
      next  # Go to the next iteration of the loop
    }
    
    # Continue with further processing if no error occurred
    permutation_category_results <- perform_overlap_permutation_analysis_multicore(
      permutations = permutations,
      seed = seed,
      num_cores = num_cores,
      reference_hgnc_vector = hgnc_symbols_vector,
      size_reference_list = row_list,
      comparison_gene_list = col_list,
      overlap_threshold = overlap_threshold,
      fdr_threshold = fdr_threshold,
      random_genes_type = random_genes_type  # Pass the random_genes_type parameter
    )
    
    permutation_category_results <- permutation_category_results %>% unlist() %>% sort()
    
    category_overlap_results$permutation_results <- permutation_category_results
    
    overlap_results[[category]] <- category_overlap_results
    
    loop_end <- Sys.time()
    loop_duration <- loop_end - loop_start
    print(loop_duration)
  }
  
  # End time measurement
  end_time <- Sys.time()
  
  # Calculate the duration
  duration <- end_time - start_time
  
  # Print the duration
  print(duration)
  
  return(overlap_results)
}
