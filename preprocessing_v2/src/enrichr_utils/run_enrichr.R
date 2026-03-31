run_enrichr <- function(gene_list, database) {
  if (!requireNamespace("enrichR", quietly = TRUE)) {
    install.packages("enrichR")
  }
  library(enrichR)
  
  # wykonanie wzbogacenia
  enriched <- enrichR::enrichr(gene_list, databases = database)
  
  # funkcja pomocnicza dodająca kolumnę n_genes
  add_n_genes <- function(df) {
    if (!"Overlap" %in% colnames(df)) return(df)
    df <- df %>%
      dplyr::mutate(
        n_genes = as.numeric(sub("/.*", "", Overlap))
      )
    return(df)
  }
  
  # jeśli baza jest pojedyncza, zwracamy df z nową kolumną
  if (length(database) == 1) {
    enriched_df <- enriched[[1]] %>% add_n_genes()
    return(enriched_df)
  } else {
    enriched <- lapply(enriched, add_n_genes)
    return(enriched)
  }
}
