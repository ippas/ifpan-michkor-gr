source("preprocessing/functions/R/install-load-packages.R")

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
