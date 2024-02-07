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

lapply(overlap_results, function(x){gene_overlap_summary(data = x)}) %>%  do.call(rbind, .) %>% 
  as.data.frame() %>% 
  # filter(permutation_FDR < 0.01) %>% 
  mutate(
    # Convert to numeric; NA if conversion fails
    number_of_significant_results = as.numeric(as.character(number_of_significant_results)),
    number_phenotype_category = as.numeric(as.character(number_phenotype_category)),
    # Replace NA with 0 after conversion
    number_of_significant_results = ifelse(is.na(number_of_significant_results), 0, number_of_significant_results),
    number_phenotype_category = ifelse(is.na(number_phenotype_category), 0, number_phenotype_category),
    # Calculate ratio; handle division by zero by replacing Inf with NA or another value
    ratio_signif_all = number_of_significant_results / number_phenotype_category
  ) %>% arrange(as.numeric(permutation_FDR)) %>% 
  filter(permutation_FDR > 0.05) %>% rownames() -> no_siginif_GR_category


lapply(categories_biobankuk, function(x){
  overlap_results[[x]]$original_data$cols
}) %>% unlist %>% unique %>% length()

lapply(categories_biobankuk, function(x){
  overlap_results[[x]]$original_data$cols
}) %>% unlist %>% unique %>% length()


read.csv("data/phenotypes/panukbiobank-phenotype-category.csv") %>%
  filter(category_type == "320_Origin_Categories") %>% .$model_name %>% unique() %>% length()

lapply(overlap_results, function(x) {
  x$significant_data$rows
}) %>% unname() %>% 
    unlist %>% length()



overlap_results$`Mental distress`$significant_data$rows

filter_phenotypes_by_category(
  genes_list = genes_phenotypes_PanUkBiobank,
  path_metafile = "data/phenotypes/panukbiobank-phenotype-category.csv",
  category_name = "Blood assays"
) %>% length()


draw_custom_heatmap(
  overlap_results$`Lifestyle and environment`,
  data_type = "significant_uniq_data",
  palette = c(
    "pastel_blue"       = "white",
    "pastel_light_blue" = "#f8dedd",
    "white"        = "#f1bcbb",
    "pastel_orange"= "#edacab",
    "pastel_red"   = "#e68a89"
  ),
  # col_mapping_vector =  clusters_mapping,
  # row_mapping_vector =  phenotyepes_mapping,
  fdr_threshold = 0.01,
  fdr_thresholds = c(0.01, 0.0001),
  color_rects =  c("#4C8D05", "#66023C"),
  color_rect = "green",
  lwd_rect = 2,
  alpha_rect = 1,
  apply_filling = F,
  color_filling = "gray",
  alpha_filling = 0.6,
  size_filling = 1,
  pch_filling = 16,
  col_significant = T,
  row_dend_width = unit(4, "cm"),  # Adjust row dendrogram width
  column_dend_height = unit(3, "cm"),  # Adjust column dendrogram height
  row_names_gp = gpar(fontsize = 12),
  column_names_gp = gpar(fontsize = 12),
  column_names_rot = 45,
  column_names_side = "top",
  overlap_threshold = 3
)


overlap_results$`Mental health`$significant_uniq_data$overlap_genes -> overlap_genes_mental_health

overlap_results$Depression$significant_uniq_data$overlap_genes -> overlap_genes_depression

overlap_results$`Psychosocial factors`$significant_uniq_data$overlap_genes -> overlap_genes_psychosocial_factors

overlap_genes_psychosocial_factors %>% length()

intersect(overlap_genes_depression, overlap_genes_psychosocial_factors) %>% length()

intersect(overlap_genes_mental_health, overlap_genes_psychosocial_factors) %>% length()

lapply(overlap_results, function(x){intersect(x$significant_uniq_data$overlap_gene, overlap_genes_psychosocial_factors) %>% length()/length(overlap_genes_psychosocial_factors)})

lapply(overlap_results, function(x){intersect(x$significant_uniq_data$overlap_gene, overlap_genes_depression) %>% length()/length(overlap_genes_depression)})

lapply(overlap_results, function(x){intersect(x$significant_uniq_data$overlap_gene, overlap_genes_mental_health) %>% length()/length(overlap_genes_mental_health)})


data.frame(
  psychosocial = sapply(overlap_results, function(x) {
    length(intersect(x$significant_uniq_data$overlap_gene, overlap_genes_psychosocial_factors)) / length(overlap_genes_psychosocial_factors)
  }),
  depression = sapply(overlap_results, function(x) {
    length(intersect(x$significant_uniq_data$overlap_gene, overlap_genes_depression)) / length(overlap_genes_depression)
  }),
  mental_health = sapply(overlap_results, function(x) {
    length(intersect(x$significant_uniq_data$overlap_gene, overlap_genes_mental_health)) / length(overlap_genes_mental_health)
  })
) -> results_df

results_df %>%
  arrange(desc(psychosocial), depression, mental_health)

