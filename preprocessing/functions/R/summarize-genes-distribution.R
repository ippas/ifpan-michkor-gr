summarize_genes_distribution <- function(data){
  Sys.time() -> start_time
  
  data$gene_lists %>% 
    unname() %>% 
    unlist() %>% 
    table %>% 
    as.data.frame() %>% 
    set_colnames(c("hgnc_symbol", "freq")) %>% 
    arrange(desc(freq)) %>% 
    mutate(total_hgnc_symbols = n()) %>%
    # filter(freq >=12) %>% 
    group_by(freq) %>%
    summarise(
      n_hgnc_symbols = n(),
      hgnc_symbol = paste(hgnc_symbol, collapse = "|"),
      proportion_hgnc_symbols = n_hgnc_symbols / first(total_hgnc_symbols)
    ) %>% 
    mutate(cumulative_sum = cumsum(proportion_hgnc_symbols)) %>% 
    ungroup %>%
    select(c(freq, n_hgnc_symbols, proportion_hgnc_symbols, cumulative_sum, hgnc_symbol)) -> genes_frequency
  
  calculate_gene_lists_diversity_metrics(data$gene_lists) -> gene_lists_diversity_metrix
  
  col_range = 1:nrow(genes_frequency)
  
  # Use sapply to iterate over the column range, ensuring output is simplified to a matrix
  intersection_number_matrix <- sapply(col_range, function(col_index) {
    # Inner sapply to iterate over each gene list and calculate intersection lengths
    sapply(data$gene_lists, function(list_item) {
      # Extract and split gene vector from the specified column and row (5th row here)
      gene_vector <- genes_frequency[col_index, "hgnc_symbol"] %>% .[[1]] %>% strsplit(., split = "\\|") %>%  unlist()
      # Calculate intersection
      intersect_genes <- intersect(list_item, gene_vector)
      # Return the length of the intersection
      length(intersect_genes)
    })
  })
  
  # Use sapply to iterate over the column range, ensuring output is simplified to a matrix
  intersection_genes_matrix <- sapply(col_range, function(col_index) {
    # Inner sapply to iterate over each gene list and calculate intersection lengths
    sapply(data$gene_lists, function(list_item) {
      # Extract and split gene vector from the specified column and row (5th row here)
      gene_vector <- genes_frequency[col_index, "hgnc_symbol"] %>% .[[1]] %>% strsplit(., split = "\\|") %>%  unlist()
      # Calculate intersection
      intersect_genes <- intersect(list_item, gene_vector)
      # Return the length of the intersection
      intersect_genes %>% paste(., collapse = "|")
    })
  })
  
  intersection_proportion_matrix <- sweep(
    intersection_number_matrix, 
    1,
    rowSums(intersection_number_matrix),
    "/") 
  
  # # calculate probability of gene in list
  # gene_propability_occurence_vector <- sapply(unname(data$gene_lists) %>% unlist() %>% unique, function(gene){
  #   calculate_gene_probability(data = data$gene_lists, gene = gene)
  # })
  
  # # calculate the gene co-occurrence matrix between pairs of genes
  # gene_co_occurrence_matrix <- calculate_gene_co_occurrence_matrix(data$gene_lists)
  # 
  # 
  # # calculate the gene probability of co-occurrence matrix between pairs of genes
  # gene_co_occurrence_probability_matrix <-
  #   calculate_gene_co_occurrence_probability_matrix(data = data$gene_lists, 
  #                                                   co_occurrence_matrix = gene_co_occurrence_matrix)
  # 
  # bayesian_probability_matrix <- sweep(gene_co_occurrence_probability_matrix, 1, gene_propability_occurence_vector, "/")
  
  # # prepare list of results
  # results <- list(
  #   genes_frequency = genes_frequency,
  #   intersection_genes_matrix = intersection_genes_matrix,
  #   intersection_number_matrix = intersection_number_matrix,
  #   intersection_proportion_matrix = intersection_proportion_matrix,
  #   gene_propability_occurence_vector = gene_propability_occurence_vector,
  #   gene_co_occurrence_matrix = gene_co_occurrence_matrix,
  #   gene_co_occurrence_probability_matrix = gene_co_occurrence_probability_matrix,
  #   bayesian_probability_matrix = bayesian_probability_matrix
  # )
  
  # prepare list of results
  results <- list(
    genes_frequency = genes_frequency,
    gene_lists_diversity_metrix = gene_lists_diversity_metrix,
    intersection_genes_matrix = intersection_genes_matrix,
    intersection_number_matrix = intersection_number_matrix,
    intersection_proportion_matrix = intersection_proportion_matrix
  )
  
  Sys.time() -> end_time
  
  end_time - start_time -> duration
  
  print(duration)
  
  return(results)
}
