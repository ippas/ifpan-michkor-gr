# Function to convert a list of gene vectors into a binary matrix
gene_list_to_matrix <- function(gene_list) {
  # Get all unique genes and sort them
  all_genes <- sort(unique(unlist(gene_list)))
  
  # Create a binary matrix with rows as genes and columns as samples
  gene_matrix <- sapply(gene_list, function(sample_genes) {
    as.integer(all_genes %in% sample_genes)
  })
  
  # Assign row names for clarity
  rownames(gene_matrix) <- all_genes
  
  # Return the binary matrix
  return(gene_matrix)
}
