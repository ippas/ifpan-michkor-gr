perform_enrichment_analysis <- function(gene_list, databases, fdr_threshold, overlap_threshold) {
  dbs_filtered <- dbs %>%
    filter(libraryName %in% databases)
  
  enriched_results <- lapply(1:length(databases), function(i){
    enrichr(
      gene_list, dbs_filtered[i,]
    )
  })
  
  names(enriched_results) <- databases
  
  enriched_results <- lapply(1:length(databases), function(i){
    enriched_results[i][[1]][3]
  }) %>% unlist(recursive = FALSE)
  
  enriched_results <- bind_rows(enriched_results, .id = "enrichr_database") %>%
    as.tibble() %>%
    mutate(n_genes = map(Genes, ~length(convert_genes_to_vector(.x, split = ";")))) %>%
    unnest(n_genes) %>%
    as.data.frame() %>%
    filter(n_genes >= overlap_threshold) %>%
    filter(Adjusted.P.value < fdr_threshold)
  
  return(enriched_results)
}
