filter_gene_list <- function(data, hgnc_symbols, drop_na = TRUE) {
  # Sprawdź czy wejście to lista
  if (!is.list(data)) {
    stop("Argument 'data' musi być listą.")
  }
  
  # Filtrowanie każdego elementu listy (wektora)
  filtered_data <- lapply(data, function(vec) {
    intersect(vec, hgnc_symbols)
  })
  
  # Opcjonalne usunięcie pustych wektorów
  if (drop_na) {
    filtered_data <- Filter(function(x) length(x) > 0, filtered_data)
  }
  
  return(filtered_data)
}
