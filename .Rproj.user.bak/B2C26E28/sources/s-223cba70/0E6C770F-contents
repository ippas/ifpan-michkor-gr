summary_overlap_results <- function(data){
  n_signif_results <- data$df %>% nrow()
  n_overlap_genes <- data$df$overlap_genes %>% convert_genes_to_vector(split = ",") %>% unique() %>% length()
  
  clusters <- data$cols
  
  lapply(c(clusters), function(cluster){
    data$df %>% 
      filter(Var2 == {{cluster}}) %>% nrow()
  }) %>% unlist -> n_signif_per_cluster
  
  lapply(c(clusters), function(cluster){
    data$df %>% 
      filter(Var2 == {{cluster}}) %>% .$overlap_genes %>% convert_genes_to_vector(split = ",") %>% unique() %>% length()
  }) %>% unlist -> n_overlap_genes_per_cluster
  
  results_per_cluster <- data.frame(
    cluster = clusters,
    n_signif_results = n_signif_per_cluster,
    n_overlap_genes = n_overlap_genes_per_cluster
  )
  
  list(
    n_signif_results = n_signif_results,
    n_overlap_genes = n_overlap_genes,
    results_per_cluster = results_per_cluster
  )
  
}

