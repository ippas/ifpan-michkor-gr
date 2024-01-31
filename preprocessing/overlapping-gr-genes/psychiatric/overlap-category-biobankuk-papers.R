source("preprocessing/functions/R/install-load-packages.R")

read_genes_from_phenotype_models(directory_path = "data/prs-models-pan-biobank-uk/") -> genes_phenotypes_PanUkBiobank



################################################################################
# Read and preprocess the data
categories_biobankuk <- read.csv("data/phenotypes/panukbiobank-phenotype-category.csv") %>%
  filter(category_type == "320_Origin_Categories") %>%
  group_by(category_name) %>%
  nest() %>%
  mutate(n = map_int(data, nrow)) %>%
  filter(n >= 10, n <= 300) %>%
  pull(category_name) 

# Initialize an empty list to store results
overlap_results <- list()

# Start time measurement
start_time <- Sys.time()

# Sequential processing using a for loop
for (category in categories_biobankuk) {
  
  print(category)
  loop_start <- Sys.time()
  
  category_genes_list <- filter_phenotypes_by_category(
    genes_list = genes_phenotypes_PanUkBiobank,
    path_metafile = "data/phenotypes/panukbiobank-phenotype-category.csv",
    category_name = category
  )
  
  # Use tryCatch to handle errors
  category_overlap_results <- tryCatch({
    analyze_gene_list_overlap(
      row_lists = papers_gene_list[!grepl("marpiech", names(papers_gene_list))],
      col_lists = category_genes_list,
      reference_hgnc_vector = hgnc_symbols_vector_v110,
      fdr_threshold = 0.01,
      overlap_threshold = 3,
      keep_original_data = FALSE
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
    permutations = 1000,
    seed = 123,
    num_cores = 30,
    reference_hgnc_vector = hgnc_symbols_vector_v110,
    size_reference_list = papers_gene_list[!grepl("marpiech", names(papers_gene_list))],
    comparison_gene_list = category_genes_list,
    overlap_threshold = 3,
    fdr_threshold = 0.01
  )
  
  permutation_category_results <- permutation_category_results %>% unlist %>% sort
  
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


################################################################################
# Set up future plan
plan(multisession, workers = 35)

# Start time measurement
start_time <- Sys.time()

# Parallel processing using future_lapply with error handling
overlap_results <- future_lapply(categories_biobankuk, function(category) {
  

  tryCatch({
    tmp <- filter_phenotypes_by_category(
      genes_list = genes_phenotypes_PanUkBiobank,
      path_metafile = "data/phenotypes/panukbiobank-phenotype-category.csv",
      category_name = category
    )
    # print(tmp)
    
    analyze_gene_list_overlap(
      row_lists = papers_gene_list[!grepl("marpiech", names(papers_gene_list))],
      col_lists = tmp,
      reference_hgnc_vector = hgnc_symbols_vector_v110,
      fdr_threshold = 0.01,
      overlap_threshold = 3,
      keep_original_data = TRUE
    )
  }, error = function(e) {
    # In case of error, return 0 or any other appropriate value
    warning(sprintf("Error in processing category '%s': %s", category, e$message))
    return(NULL)
  })

}, future.seed = TRUE)

names(overlap_results) <- categories_biobankuk

# End time measurement
end_time <- Sys.time()

# Calculate the duration
duration <- end_time - start_time

# Print the duration
print(duration)

Filter(Negate(is.null), overlap_results) -> overlap_results
overlap_results$`Mental distress`$original_data$df %>% 
  filter(gene_overlap_count > 1) %>% 
  filter(fdr < 0.9)
  mutate(log10_fdr = log10(fdr))


lapply(overlap_results, function(x) {x$significant_uniq_data$df %>% nrow})

overlap_results$Anxiety$significant_uniq_data$df
