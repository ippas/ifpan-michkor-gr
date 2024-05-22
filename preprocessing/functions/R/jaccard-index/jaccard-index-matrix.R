# Function to calculate Jaccard index between two sets
jaccard_index <- function(set1, set2) {
  if (length(set1) == 0 & length(set2) == 0) {
    return(1)  # if both sets are empty, define their similarity as 1
  }
  intersection_size <- length(intersect(set1, set2))
  union_size <- length(union(set1, set2))
  return(intersection_size / union_size)
}

# Function to create a Jaccard index matrix from a list of gene vectors
jaccard_matrix <- function(gene_list) {
  num_genes <- length(gene_list)
  matrix_result <- matrix(0, nrow = num_genes, ncol = num_genes)
  for (i in 1:num_genes) {
    for (j in i:num_genes) {
      if (i == j) {
        matrix_result[i, j] <-NA
      } else {
        index <- jaccard_index(gene_list[[i]], gene_list[[j]])
        matrix_result[i, j] <- index
        matrix_result[j, i] <- index  # since the matrix is symmetric
      }
    }
  }
  colnames(matrix_result) <- rownames(matrix_result) <- names(gene_list)
  return(matrix_result)
}


jaccard_index_per_tissue <- function(gene_list, tissue_data, regulation, size_list, rank_criterion){
  tissue_data %>% 
    filter(regulation == {{regulation}},
           size_list == {{size_list}},
           rank_criterion == {{rank_criterion}}) -> preprocessing_data
  
  preprocessing_data$tissue %>% unique -> tissue_vector
  
  lapply(tissue_vector, function(tissue_name) {
    jaccard_index(set1 = gene_list, 
                  set2 = preprocessing_data %>% filter(tissue == tissue_name) %>% .$hgnc_symbol)
  }) -> jaccard_results
  
  preprocessing_data %>% 
    select(c(tissue, n_lists_filt)) %>% 
    unique() %>% 
    mutate(name = paste(tissue, n_lists_filt, sep = " ")) %>% .$name %>% unique -> names_column
  
  print(names_column)
   
  names(jaccard_results) <- names_column
  
  jaccard_results %>% as.data.frame() -> jaccard_results
  
  return(jaccard_results)
}

regulation_norm_master_df


jaccard_index_multiple_comparison <- function(regulation_master_data,
                                              tissue_data,
                                              parameters_df){
  tissue_data_2 <- tissue_data
  
  parameters_df %>% 
    mutate(name = paste(regulation, size_list, rank_criterion, paste0("nListPlus", additional_adjustment), sep = "_")) -> parameters_df
  
  
  lapply(c(1:nrow({{parameters_df}})), function(i){
    print(parameters_df[i,])
    
    gene_list_2 <- regulation_master_data[[as.character(parameters_df[i, "name"])]]$gene_list
    print(gene_list_2)
    
    
    jaccard_index_per_tissue(gene_list = gene_list_2,
                             tissue_data = tissue_data_2,
                             regulation = parameters_df[i, "regulation"] %>% as.character(),
                             size_list = parameters_df[i, "size_list"] %>% as.numeric,
                             rank_criterion = parameters_df[i, "rank_criterion"] %>% as.character()) -> result
    result
    
  }) -> jaccard_results
  
  
  names(jaccard_results) <- parameters_df$name
  
  return(jaccard_results) 
}
