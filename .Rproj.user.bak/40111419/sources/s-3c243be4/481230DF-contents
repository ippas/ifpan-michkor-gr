analyze_gene_list_overlap <-
  function(row_lists,
           col_lists,
           reference_hgnc_vector,
           fdr_threshold = 0.05,
           overlap_threshold = 3,
           keep_original_data = TRUE) {
    # Combine and get unique elements from both lists
    combined_list <- lapply(c(row_lists, col_lists), unique)
    
    # Perform chi2 tests or similar analysis
    analysis_results <-
      perform_chi2_tests(combined_list, reference_hgnc_vector)
    
    # Process overlap results or similar post-analysis processing
    processed_data <- processing_overlap_results(
      data = analysis_results,
      rows_to_filter = !rownames(analysis_results$p_value_matrix) %in% names(row_lists),
      cols_to_filter = names(row_lists),
      overlap_threshold = {{overlap_threshold}},
      fdr_threshold = {{fdr_threshold}},
      genes_list = combined_list
    )
    
    
    # Remove the "original_data" element from the list
    if(!keep_original_data){
      processed_data <- processed_data[!(names(processed_data) == "original_data")]
    }
    
    return(processed_data)
  }
