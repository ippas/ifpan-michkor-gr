gene_overlap_summary <- function(data) {
  # Calculate the number of significant results
  signif_overlap <- nrow(data$significant_uniq_data$df)
  
  signif_phenotypes <- data$significant_uniq_data$df$Var1 %>%  unique() %>% length()
  signif_papers <- data$significant_uniq_data$df$Var2 %>%  unique() %>% length()
  # print(number_of_significant_results)
  
  # Calculate permutation FDR
  permutation_FDR <- sum(data$permutation_results >= signif_overlap) / 1000
  # print(permutation_FDR)
  
  max_permutation <- max(data$permutation_results)
  
  df_genes_freq <-
    data$significant_uniq_data$df$overlap_genes %>% lapply(., function(x) {
      strsplit(x, split = ",") %>% unlist
    }) %>% unlist %>% table %>% as.data.frame() %>% set_colnames(c("gene_name", "freq")) %>%  arrange(desc(freq)) %>% 
    filter(freq > 1)
  
  data$original_data$cols %>% length() -> number_phenotype_category
  
  # Create a summary object
  summary <- list(
    signif_overlap = signif_overlap,
    signif_phenotypes = signif_phenotypes,
    signif_papers = signif_papers,
    permutation_FDR = permutation_FDR,
    max_permutation_signif_overlap = max_permutation,
    number_phenotype_category = number_phenotype_category
  )
  
  # Return the summary object
  return(summary)
}
