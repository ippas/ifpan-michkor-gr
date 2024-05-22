calculate_gene_probability <- function(data, gene) {
  # Calculate the total number of cases
  total_cases <- length(data)
  
  # Count the number of cases that include the specified gene
  cases_with_gene <- sum(sapply(data, function(case) gene %in% case))
  
  # Calculate the probability of the gene being involved
  probability_gene <- cases_with_gene / total_cases
  
  # Return the probability
  return(probability_gene)
}

# Define the function to calculate the co-occurrence matrix
calculate_gene_co_occurrence_matrix <- function(gene_data) {
  Sys.time() -> start_time
  # Flatten the list to create edges
  edges <- unlist(lapply(gene_data, function(set) {
    combn(set, 2, simplify = FALSE)
  }), recursive = FALSE)
  
  # Unlist the pairs to a single vector
  edges <- unlist(edges)
  
  # Create a graph from the edge list
  g <- graph(edges, directed = FALSE)
  
  # Generate the adjacency (co-occurrence) matrix
  co_occurrence_matrix <- as_adjacency_matrix(g, type = "both", attr = NULL, sparse = FALSE)
  
  Sys.time() -> end_time
  
  print(end_time - start_time)
  
  return(co_occurrence_matrix)
}


# Function to calculate gene co-occurrence probability matrix
calculate_gene_co_occurrence_probability_matrix <- function(data, co_occurrence_matrix){
  # Args:
  #   data: A list or collection of gene sets
  #   co_occurrence_matrix: A matrix representing the co-occurrence counts between pairs of genes
  # Returns:
  #   A matrix representing the probability of co-occurrence between pairs of genes
  
  # Calculate the total number of gene lists (sets)
  number_of_list <- length(data)
  
  # Normalize the co-occurrence matrix by the number of lists to get probabilities
  # This converts raw counts to probabilities by dividing by the total number of lists
  co_occurrence_probability_matrix <- co_occurrence_matrix / number_of_list
  
  # Set the diagonal elements to 0 to ignore self-co-occurrences
  # This is common when focusing on the relationship between different genes
  diag(co_occurrence_probability_matrix) <- 0
  
  # Return the normalized co-occurrence probability matrix
  return(co_occurrence_probability_matrix)
}
