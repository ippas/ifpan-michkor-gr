generate_random_gene_lists <- function(gene_list, n_lists, length) {
  if (length > length(gene_list)) {
    stop("❌ length cannot exceed the number of available genes in gene_list")
  }
  
  random_sets <- vector("list", n_lists)
  
  for (i in seq_len(n_lists)) {
    random_sets[[i]] <- sample(gene_list, size = length, replace = FALSE)
  }
  
  names(random_sets) <- paste0("randomGeneList_", seq_len(n_lists))
  return(random_sets)
}
