refine_gene_lists <- function(data, size_threshold = 0, columns, cumsum_thresholds, freq_thresholds, keep_column = NULL) {
  
  # Process the input data
  refined_lists <- refine_and_label_gene_lists(data, columns, size_threshold = size_threshold, keep_column = keep_column)
  genes_distribution <- summarize_genes_distribution(refined_lists)
  master_gene_lists <- create_master_gene_lists(cumsum_thresholds, genes_distribution$genes_frequency)
  rare_gene_lists <- create_rare_gene_lists(freq_thresholds, genes_distribution$genes_frequency)
  
  # Return the processed lists
  list(
    refine_gene_lists = refined_lists,
    genes_distribution = genes_distribution,
    master_gene_lists = master_gene_lists,
    rare_gene_lists = rare_gene_lists
  )
}
