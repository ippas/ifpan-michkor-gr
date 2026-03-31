processing_overlap_results_fisher <- function(data, genes_list, rows_to_filter, cols_to_filter,
                                              fdr_threshold = 0.05, overlap_threshold = 2) {
  message("Start: processing_overlap_results_fisher()")
  
  # Step 1: Extract matrices
  message("Extracting and reshaping input matrices...")
  pval_mat     <- data$p_value_matrix[rows_to_filter, cols_to_filter, drop = FALSE]
  fisher_mat   <- data$fisher_value_matrix[rows_to_filter, cols_to_filter, drop = FALSE]
  overlap_mat  <- data$overlap_genes_matrix[rows_to_filter, cols_to_filter, drop = FALSE]
  count_mat    <- data$number_overlap_matrix[rows_to_filter, cols_to_filter, drop = FALSE]
  
  # Melt each matrix
  df_pval     <- as.data.frame(as.table(pval_mat), stringsAsFactors = FALSE)
  df_fisher   <- as.data.frame(as.table(fisher_mat), stringsAsFactors = FALSE)
  df_overlap  <- as.data.frame(as.table(overlap_mat), stringsAsFactors = FALSE)
  df_count    <- as.data.frame(as.table(count_mat), stringsAsFactors = FALSE)
  
  colnames(df_pval)     <- c("Var1", "Var2", "p_value")
  colnames(df_fisher)   <- c("Var1", "Var2", "fisher_value")
  colnames(df_overlap)  <- c("Var1", "Var2", "overlap_genes")
  colnames(df_count)    <- c("Var1", "Var2", "gene_overlap_count")
  
  # Step 2: Merge and process
  original_df <- df_pval %>%
    dplyr::left_join(df_fisher,   by = c("Var1", "Var2")) %>%
    dplyr::left_join(df_overlap,  by = c("Var1", "Var2")) %>%
    dplyr::left_join(df_count,    by = c("Var1", "Var2")) %>%
    dplyr::mutate(
      Var1 = as.character(Var1),
      Var2 = as.character(Var2),
      overlap_genes = ifelse(is.na(overlap_genes), "", overlap_genes),
      gene_overlap_count = ifelse(is.na(gene_overlap_count), 0, gene_overlap_count),
      fdr = p.adjust(p_value, method = "fdr"),
      overlap_genes = strsplit(overlap_genes, ",") %>%
        purrr::map(~sort(.) %>% paste(collapse = ",")) %>%
        unlist()
    )
  
  message(paste("Extracted", nrow(original_df), "rows"))
  
  # Step 3: FDR matrix
  fdr_value_matrix <- reshape2::melt(pval_mat) %>%
    dplyr::mutate(value = p.adjust(value, method = "fdr")) %>%
    reshape2::dcast(Var1 ~ Var2, value.var = "value") %>%
    tibble::column_to_rownames("Var1") %>%
    as.matrix()
  
  original_list <- list(
    p_value_matrix = pval_mat,
    fisher_value_matrix = fisher_mat,
    overlap_genes_matrix = overlap_mat,
    number_overlap_matrix = count_mat,
    fdr_value_matrix = fdr_value_matrix
  )
  
  original_rows <- unique(original_df$Var1)
  original_cols <- unique(original_df$Var2)
  original_overlap_genes <- unique(unlist(strsplit(original_df$overlap_genes, ",")))
  
  # Step 4: Significant
  message("Filtering significant results...")
  significant_df <- original_df %>%
    dplyr::filter(fdr < fdr_threshold, gene_overlap_count >= overlap_threshold)
  
  if (nrow(significant_df) == 0) {
    warning("No significant results found.")
    significant_df <- original_df[0, ]
  }
  
  significant_rows <- unique(significant_df$Var1)
  significant_cols <- unique(significant_df$Var2)
  significant_overlap_genes <- unique(unlist(strsplit(significant_df$overlap_genes, ",")))
  
  significant_list <- list(
    p_value_matrix = pval_mat[significant_rows, cols_to_filter, drop = FALSE],
    fisher_value_matrix = fisher_mat[significant_rows, cols_to_filter, drop = FALSE],
    overlap_genes_matrix = overlap_mat[significant_rows, cols_to_filter, drop = FALSE],
    number_overlap_matrix = count_mat[significant_rows, cols_to_filter, drop = FALSE],
    fdr_value_matrix = fdr_value_matrix[significant_rows, cols_to_filter, drop = FALSE]
  )
  
  # Step 5: Unique
  message("Filtering uniquely significant results...")
  significant_uniq_df <- significant_df %>%
    dplyr::group_by(overlap_genes, Var2) %>%
    dplyr::slice_min(fdr, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(Var1, Var2, p_value, fisher_value, gene_overlap_count, overlap_genes, fdr)
  
  significant_uniq_rows <- unique(significant_uniq_df$Var1)
  significant_uniq_cols <- unique(significant_uniq_df$Var2)
  significant_uniq_overlap_genes <- unique(unlist(strsplit(significant_uniq_df$overlap_genes, ",")))
  
  significant_uniq_list <- list(
    p_value_matrix = pval_mat[significant_uniq_rows, cols_to_filter, drop = FALSE],
    fisher_value_matrix = fisher_mat[significant_uniq_rows, cols_to_filter, drop = FALSE],
    overlap_genes_matrix = overlap_mat[significant_uniq_rows, cols_to_filter, drop = FALSE],
    number_overlap_matrix = count_mat[significant_uniq_rows, cols_to_filter, drop = FALSE],
    fdr_value_matrix = fdr_value_matrix[significant_uniq_rows, cols_to_filter, drop = FALSE]
  )
  
  # Step 6: Gene list sizes
  gene_list_sizes <- sapply(genes_list, length)
  
  # Step 7: Return
  return(list(
    original_data = list(
      list = original_list,
      df = original_df,
      rows = original_rows,
      cols = original_cols,
      overlap_genes = original_overlap_genes
    ),
    significant_data = list(
      list = significant_list,
      df = significant_df,
      rows = significant_rows,
      cols = significant_cols,
      overlap_genes = significant_overlap_genes
    ),
    significant_uniq_data = list(
      list = significant_uniq_list,
      df = significant_uniq_df,
      rows = significant_uniq_rows,
      cols = significant_uniq_cols,
      overlap_genes = significant_uniq_overlap_genes
    ),
    gene_list_sizes = gene_list_sizes
  ))
}
