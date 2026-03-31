source("preprocessing/functions/R/install-load-packages.R")


read_genes_from_phenotype_models(directory_path = "data/prs-models-pan-biobank-uk/") -> genes_phenotypes_PanUkBiobank

filter_phenotypes_by_category(
  genes_list = genes_phenotypes_PanUkBiobank,
  path_metafile = "data/phenotypes/panukbiobank-phenotype-category.csv",
  category_name = "Mental health"
) -> mental_health_genes_list


Sys.time() -> start_time
permutation_mental_health_results2 <- perform_overlap_permutation_analysis_parallel(
  permutations = 100,  # Number of permutations
  seed = 123,          # Initial seed for reproducibility
  num_cores = 35,       # Number of cores to use
  # Additional arguments for analyze_random_gene_sets
  reference_hgnc_vector = hgnc_symbols_vector_v110,
  size_reference_list = papers_gene_list[!grepl("marpiech", names(papers_gene_list))],
  comparison_gene_list = mental_health_genes_list,
  overlap_threshold = 3,
  fdr_threshold = 0.01
)

Sys.time() -> end_time

start_time - end_time


# Now, use the save function to save it to an .RData file
save(permutation_mental_health_results, file = "results/permutation-mental-health-results.RData")

load("results/permutation-mental-health-results.RData")

permutation_mental_health_results

################################################################################
filter_phenotypes_by_category(
  genes_list = genes_phenotypes_PanUkBiobank,
  path_metafile = "data/phenotypes/panukbiobank-phenotype-category.csv",
  category_name = "Depression"
) -> depression_category_genes_list

Sys.time() -> start_time
permutation_depression_results <- perform_overlap_permutation_analysis_parallel(
  permutations = 1000,  # Number of permutations
  seed = 123,          # Initial seed for reproducibility
  num_cores = 35,       # Number of cores to use
  # Additional arguments for analyze_random_gene_sets
  reference_hgnc_vector = hgnc_symbols_vector_v110,
  size_reference_list = papers_gene_list[!grepl("marpiech", names(papers_gene_list))],
  comparison_gene_list = depression_category_genes_list,
  overlap_threshold = 3,
  fdr_threshold = 0.01
)

Sys.time() -> end_time

start_time - end_time

# Now, use the save function to save it to an .RData file
save(permutation_depression_results, file = "results/permutation_depression_results.RData")


################################################################################
filter_phenotypes_by_category(
  genes_list = genes_phenotypes_PanUkBiobank,
  path_metafile = "data/phenotypes/panukbiobank-phenotype-category.csv",
  category_name = "Psychosocial factors"
) ->  psychosocial_factors_genes_list

Sys.time() -> start_time
permutation_psychosocial_factors_results <- perform_overlap_permutation_analysis_parallel(
  permutations = 10,  # Number of permutations
  seed = 123,          # Initial seed for reproducibility
  num_cores = 35,       # Number of cores to use
  # Additional arguments for analyze_random_gene_sets
  reference_hgnc_vector = hgnc_symbols_vector_v110,
  size_reference_list = papers_gene_list[!grepl("marpiech", names(papers_gene_list))],
  comparison_gene_list = psychosocial_factors_genes_list,
  overlap_threshold = 3,
  fdr_threshold = 0.01
)

Sys.time() -> end_time

start_time - end_time

# Now, use the save function to save it to an .RData file
save(permutation_psychosocial_factors_results, file = "results/permutation_psychosocial_factors_results.RData")
