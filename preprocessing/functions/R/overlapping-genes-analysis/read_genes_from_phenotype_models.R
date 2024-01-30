read_genes_from_phenotype_models <- function(directory_path, depression_models_pattern) {
  list.files(directory_path, pattern = "\\.ya?ml$", full.names = TRUE) %>%
    # keep(~ grepl(depression_models_pattern, .x)) %>% 
    set_names(map(., ~ basename(.x) %>% tools::file_path_sans_ext())) %>%
    # Read YAML content from each file
    map(read_yaml) %>%
    # Extract genes information and name the list elements after the files
    map(~ .x$description$genes) %>% 
    lapply(., unique)
}
