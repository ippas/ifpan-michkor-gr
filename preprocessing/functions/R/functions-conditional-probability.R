# Define a function to calculate the probability of a given gene being involved
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

# 
# Wrapper function to calculate the matrix of conditional probabilities
calculate_conditional_probability_matrix <- function(data) {
  # This function calculates the conditional probability of observing gene2 given gene1
  # within a list of gene lists.
  #
  # Args:
  #   data: A list of character vectors, where each vector contains gene names
  #   gene1: The gene condition (the gene of interest whose presence is given)
  #   gene2: The outcome gene (the gene whose conditional probability is being calculated)
  #
  # Returns:
  #   The conditional probability of gene2 given gene1.
  # Extract unique genes from the data
  unique_genes <- unique(unlist(data))
  
  # Initialize an empty matrix to store the probabilities
  probability_matrix <- matrix(0, nrow = length(unique_genes), ncol = length(unique_genes),
                               dimnames = list(unique_genes, unique_genes))
  
  # Populate the matrix with conditional probabilities
  for (i in 1:length(unique_genes)) {
    for (j in 1:length(unique_genes)) {
      # if (i != j) { # To avoid calculating probability of a gene given itself, unless desired
      probability_matrix[i, j] <- calculate_conditional_probability(data, unique_genes[i], unique_genes[j])
      # }
    }
  }
  
  attr(probability_matrix, "descriptions") <- c(gene_in_rows = "gene1: The gene condition (the gene of interest whose presence is given)",
                                                gene_in_columns = "gene2: The outcome gene (the gene whose conditional probability is being calculated)")
  
  return(probability_matrix)
}
