filter_phenotypes_by_category <- function(genes_list, path_metafile, category_name) {
  # Read the metafile
  phenotype_models <- read.csv(path_metafile) %>%
    filter(category_name == {{category_name}}) %>%
    .$model_name %>%
    unique()
  
  # Use lapply to filter genes list for each phenotype model
  filtered_genes_list <- lapply(phenotype_models, function(model) {
    genes_list[grepl(model, names(genes_list))]
  }) %>% unlist(recursive = F)
  
  return(filtered_genes_list)
}

  