generate_random_gene_lists_by_lengths <- function(gene_list, lengths) {
  
  # --- sanity checks ---
  if (!is.numeric(lengths) || any(lengths <= 0)) {
    stop("❌ lengths must be a numeric vector of positive integers")
  }
  
  if (any(lengths > length(gene_list))) {
    stop("❌ one or more requested lengths exceed the number of available genes")
  }
  
  lengths <- as.integer(lengths)
  
  # --- generate random gene sets ---
  random_sets <- lapply(lengths, function(len) {
    sample(gene_list, size = len, replace = FALSE)
  })
  
  names(random_sets) <- paste0(
    "randomGeneList_", 
    seq_along(lengths), 
    "_n", lengths
  )
  
  return(random_sets)
}
