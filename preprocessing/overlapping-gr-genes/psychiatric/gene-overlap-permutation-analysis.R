source("preprocessing/functions/R/install-load-packages.R")

source("preprocessing/functions/R/gene-overlap-permutation-analysis/generate-random-gene-sets.R")
source("preprocessing/functions/R/gene-overlap-permutation-analysis/analyze-random-gene-sets.R")
source("preprocessing/functions/R/gene-overlap-permutation-analysis/perform-overlap-permutation-analysis-parallel.R")
source("preprocessing/functions/R/perform-chi2-tests.R")


read.csv("data/phenotypes/panukbiobank-phenotype-category.csv") %>% 
  filter(category_name == "Mental health") %>% 
  .$model_name %>% unique %>%  paste(., collapse = "|") -> depression_models_pattern

# List all files in the directory with full names
list.files("data/prs-models-pan-biobank-uk/", full.names = TRUE) -> model_files_vector

model_files_vector %>% 
  keep(~ grepl(depression_models_pattern, .x)) %>% 
  set_names(map(., ~ basename(.x) %>% tools::file_path_sans_ext())) %>%
  # Read YAML content from each file
  map(read_yaml) %>%
  # Extract genes information and name the list elements after the files
  map(~ .x$description$genes) %>% 
  lapply(., unique) -> depression_genes_list

Sys.time() -> start_time
permutation_mental_health_results2 <- perform_overlap_permutation_analysis_parallel(
  permutations = 10,  # Number of permutations
  seed = 223,          # Initial seed for reproducibility
  num_cores = 35,       # Number of cores to use
  # Additional arguments for analyze_random_gene_sets
  reference_hgnc_vector = hgnc_symbols_vector_v110,
  size_reference_list = papers_gene_list[!grepl("marpiech", names(papers_gene_list))],
  comparison_gene_list = depression_genes_list,
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
read.csv("data/phenotypes/panukbiobank-phenotype-category.csv") %>% 
  filter(category_name == "Depression") %>% 
  .$model_name %>% unique %>%  paste(., collapse = "|") -> depression_models_pattern

# List all files in the directory with full names
list.files("data/prs-models-pan-biobank-uk/", full.names = TRUE) -> model_files_vector

model_files_vector %>% 
  keep(~ grepl(depression_models_pattern, .x)) %>% 
  set_names(map(., ~ basename(.x) %>% tools::file_path_sans_ext())) %>%
  # Read YAML content from each file
  map(read_yaml) %>%
  # Extract genes information and name the list elements after the files
  map(~ .x$description$genes) %>% 
  lapply(., unique) -> depression_category_genes_list

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
read.csv("data/phenotypes/panukbiobank-phenotype-category.csv") %>% 
  filter(category_name == "Psychosocial factors") %>% 
  .$model_name %>% unique %>%  paste(., collapse = "|") -> psychosocial_factors_models_pattern


# List all files in the directory with full names
list.files("data/prs-models-pan-biobank-uk/", full.names = TRUE) -> model_files_vector

model_files_vector %>% 
  keep(~ grepl(psychosocial_factors_models_pattern, .x)) %>% 
  set_names(map(., ~ basename(.x) %>% tools::file_path_sans_ext())) %>%
  # Read YAML content from each file
  map(read_yaml) %>%
  # Extract genes information and name the list elements after the files
  map(~ .x$description$genes) %>% 
  lapply(., unique) -> psychosocial_factors_genes_list


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
