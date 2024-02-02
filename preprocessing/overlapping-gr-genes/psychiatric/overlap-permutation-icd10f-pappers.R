source("preprocessing/functions/R/install-load-packages.R")

filter_phenotypes_by_category(
  genes_list = genes_phenotypes_PanUkBiobank,
  path_metafile = "data/phenotypes/panukbiobank-phenotype-category.csv",
  category_name = "Psychosocial factors"
) ->  psychosocial_factors_genes_list

read.csv("data/phenotypes/panukbiobank-phenotype-category.csv") %>% 
  filter(field_id == "41270") %>% 
  filter(category_type == "320_Origin_Categories") %>%
  select(model_name) %>% 
  unique() %>% 
  mutate(icd10 = sapply(str_split(model_name, "-"), `[`, 2)) %>% 
  filter(grepl("^[A-Za-z]", icd10)) %>% 
  mutate(category_name = substr(icd10, 1, 1)) %>% 
  group_by(category_name) %>%
  mutate(count_icd10_category = n()) %>%
  ungroup() %>% 
  write.csv(., file = "data/phenotypes/panukbiobank-phenotype-icd10-category.csv", col.names = TRUE, row.names = FALSE, quote = FALSE)

read.csv(file = "data/phenotypes/panukbiobank-phenotype-icd10-category.csv") %>% 
  filter(count_icd10_category > 10) %>% .$category_name %>% unique -> icd10_categories_biobankuk

# Initialize an empty list to store results
overlap_icd10_category_results <- list()

# Start time measurement
start_time <- Sys.time()

# Sequential processing using a for loop
for (category in icd10_categories_biobankuk) {
  
  print(category)
  loop_start <- Sys.time()
  
  category_genes_list <- filter_phenotypes_by_category(
    genes_list = genes_phenotypes_PanUkBiobank,
    path_metafile = "data/phenotypes/panukbiobank-phenotype-icd10-category.csv",
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
      keep_original_data = TRUE
    )
  }, error = function(e) {
    warning(sprintf("Error in analyze_gene_list_overlap for category '%s': %s", category, e$message))
    return(NULL)  # Return NULL in case of error
  })
  
  # Skip further processing ifoverlap_icd10_category_resultsis NULL
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
  
  overlap_icd10_category_results[[category]] <- category_overlap_results
  
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

lapply(overlap_icd10_category_results, function(x){gene_overlap_summary(data = x)}) %>%  do.call(rbind, .) %>% 
  as.data.frame() %>% 
  filter(permutation_FDR < 0.1) %>% 
  mutate(
    # Convert to numeric; NA if conversion fails
    number_of_significant_results = as.numeric(as.character(number_of_significant_results)),
    number_phenotype_category = as.numeric(as.character(number_phenotype_category)),
    # Replace NA with 0 after conversion
    number_of_significant_results = ifelse(is.na(number_of_significant_results), 0, number_of_significant_results),
    number_phenotype_category = ifelse(is.na(number_phenotype_category), 0, number_phenotype_category),
    # Calculate ratio; handle division by zero by replacing Inf with NA or another value
    ratio_signif_all = number_of_significant_results / number_phenotype_category
  ) 


overlap_icd10_category_results$c$original_data$cols %>% length()













read.csv("data/phenotypes/panukbiobank-phenotype-category.csv") %>% 
  filter(field_id == "41270") %>% 
  filter(category_type == "320_Origin_Categories") %>%
  select(model_name) %>% 
  unique() %>% 
  mutate(icd10 = sapply(str_split(model_name, "-"), `[`, 2)) %>% 
  filter(grepl("f", icd10)) %>% 
  .$model_name %>% unique() %>% paste(., collapse = "|") -> icd10f_models_pattern

# List all files in the directory with full names
list.files("data/prs-models-pan-biobank-uk/", full.names = TRUE) -> model_files_vector

model_files_vector %>% 
  keep(~ grepl(icd10f_models_pattern, .x)) %>% 
  set_names(map(., ~ basename(.x) %>% tools::file_path_sans_ext())) %>%
  # Read YAML content from each file
  map(read_yaml) %>%
  # Extract genes information and name the list elements after the files
  map(~ .x$description$genes) %>% 
  lapply(., unique) -> icd10f_genes_list


Sys.time() -> start_time
permutation_icd10f_results <- perform_overlap_permutation_analysis_parallel(
  permutations = 1000,  # Number of permutations
  seed = 123,          # Initial seed for reproducibility
  num_cores = 35,       # Number of cores to use
  # Additional arguments for analyze_random_gene_sets
  reference_hgnc_vector = hgnc_symbols_vector_v110,
  size_reference_list = papers_gene_list[!grepl("marpiech", names(papers_gene_list))],
  comparison_gene_list = icd10f_genes_list,
  overlap_threshold = 3,
  fdr_threshold = 0.01
)

Sys.time() -> end_time

start_time - end_time


  
