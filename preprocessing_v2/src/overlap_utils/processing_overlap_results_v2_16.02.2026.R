processing_overlap_results <- function(
    data,
    genes_list,
    rows_to_filter,
    cols_to_filter,
    fdr_threshold = 0.05,
    overlap_threshold = 2,
    total_genes = NULL,     # potrzebne do OR i expected (fallback)
    epsilon = 1e-6          # pseudo-count do OR i do FE (jeśli kiedyś dodasz)
) {
  # ============================================================
  # processing_overlap_results (v2)
  # - propagates expected_overlap_matrix into df (expected_overlap)
  # - adds Var1_n_genes and Var2_n_genes to df
  # - computes OR/log2OR if missing (as before)
  # - computes expected_overlap_matrix if missing (fallback)
  # ============================================================
  
  # ---------- packages ----------
  library(dplyr)
  library(tidyr)
  library(reshape2)
  library(purrr)
  library(tibble)
  
  # ---------- (0) basic checks ----------
  if (is.null(data$number_overlap_matrix)) {
    stop("`data$number_overlap_matrix` is required.")
  }
  overlap_mat <- data$number_overlap_matrix
  
  if (is.null(rownames(overlap_mat)) || is.null(colnames(overlap_mat))) {
    stop("`number_overlap_matrix` must have rownames and colnames.")
  }
  
  rn_all <- rownames(overlap_mat)
  cn_all <- colnames(overlap_mat)
  
  # check that genes_list is named
  if (is.null(names(genes_list)) || any(names(genes_list) == "")) {
    stop("`genes_list` must be a named list. Names must match matrix dimnames.")
  }
  
  # sizes for Var1/Var2
  gene_list_sizes <- sapply(genes_list, length)
  if (!all(rn_all %in% names(gene_list_sizes))) {
    missing <- setdiff(rn_all, names(gene_list_sizes))
    stop("Missing gene lists for matrix rownames in `genes_list`: ", paste(missing, collapse = ", "))
  }
  if (!all(cn_all %in% names(gene_list_sizes))) {
    missing <- setdiff(cn_all, names(gene_list_sizes))
    stop("Missing gene lists for matrix colnames in `genes_list`: ", paste(missing, collapse = ", "))
  }
  
  # ---------- (A) If OR/log2OR are missing, compute them ----------
  if (is.null(data$odds_ratio_matrix) || is.null(data$log2_odds_ratio_matrix)) {
    if (is.null(total_genes)) {
      stop("`total_genes` is required to compute odds ratios when they are missing in `data`.")
    }
    
    n_bg <- length(total_genes)
    
    rn <- rownames(overlap_mat); cn <- colnames(overlap_mat)
    or_mat <- matrix(NA_real_, nrow = nrow(overlap_mat), ncol = ncol(overlap_mat),
                     dimnames = list(rn, cn))
    log2or_mat <- or_mat
    
    for (i in seq_along(rn)) {
      for (j in seq_along(cn)) {
        a <- overlap_mat[i, j]
        ni <- gene_list_sizes[rn[i]]
        nj <- gene_list_sizes[cn[j]]
        
        # standard 2x2:
        # a = overlap
        # b = in i only
        # c = in j only
        # d = in neither (background)
        b <- ni - a
        c <- nj - a
        d <- n_bg - (ni + nj - a)
        
        or <- ((a + epsilon) * (d + epsilon)) / ((b + epsilon) * (c + epsilon))
        or_mat[i, j] <- or
        log2or_mat[i, j] <- log2(or)
      }
    }
    
    data$odds_ratio_matrix <- or_mat
    data$log2_odds_ratio_matrix <- log2or_mat
  }
  
  # ---------- (A2) If expected_overlap_matrix is missing, compute it ----------
  if (is.null(data$expected_overlap_matrix)) {
    if (is.null(total_genes)) {
      stop("`total_genes` is required to compute expected overlaps when they are missing in `data`.")
    }
    
    n_bg <- length(total_genes)
    rn <- rownames(overlap_mat); cn <- colnames(overlap_mat)
    
    exp_mat <- matrix(NA_real_, nrow = nrow(overlap_mat), ncol = ncol(overlap_mat),
                      dimnames = list(rn, cn))
    
    for (i in seq_along(rn)) {
      for (j in seq_along(cn)) {
        ni <- gene_list_sizes[rn[i]]
        nj <- gene_list_sizes[cn[j]]
        exp_mat[i, j] <- (ni * nj) / n_bg
      }
    }
    
    data$expected_overlap_matrix <- exp_mat
  }
  
  # ---------- (B) internal helpers ----------
  safe_melt <- function(x, nm) {
    if (is.null(x)) return(NULL)
    x[rows_to_filter, cols_to_filter, drop = FALSE] %>%
      melt() %>%
      `colnames<-`(c("Var1", "Var2", nm))
  }
  
  filter_matrices <- function(results) {
    out <- list()
    keep <- intersect(
      names(results),
      c(
        "p_value_matrix",
        "chi2_value_matrix",
        "number_overlap_matrix",
        "expected_overlap_matrix",      # NEW ✅
        "overlap_genes_matrix",
        "odds_ratio_matrix",
        "log2_odds_ratio_matrix",
        "fdr_matrix"
      )
    )
    for (nm in keep) {
      out[[nm]] <- results[[nm]][rows_to_filter, cols_to_filter, drop = FALSE]
    }
    out
  }
  
  # ---------- (C) Build original_df ----------
  pieces <- list(
    safe_melt(data$p_value_matrix,          "p_value"),
    safe_melt(data$chi2_value_matrix,       "chi2"),
    safe_melt(data$number_overlap_matrix,   "gene_overlap_count"),
    safe_melt(data$expected_overlap_matrix, "expected_overlap"),   # NEW ✅
    safe_melt(data$overlap_genes_matrix,    "overlap_genes"),
    safe_melt(data$odds_ratio_matrix,       "odds_ratio"),
    safe_melt(data$log2_odds_ratio_matrix,  "log2_odds_ratio"),
    safe_melt(data$fdr_matrix,              "fdr_value")
  )
  pieces <- pieces[!vapply(pieces, is.null, logical(1))]
  
  original_df <- Reduce(
    function(x, y) dplyr::left_join(x, y, by = c("Var1", "Var2")),
    pieces
  )
  
  # ---------- (C2) Add sizes for Var1/Var2 + FDR + normalize overlap_genes ----------
  original_df <- original_df %>%
    mutate(
      Var1 = as.character(Var1),
      Var2 = as.character(Var2),
      
      Var1_n_genes = unname(gene_list_sizes[Var1]),  # NEW ✅
      Var2_n_genes = unname(gene_list_sizes[Var2]),  # NEW ✅
      
      # global FDR across all melted tests (as in your original)
      fdr = p.adjust(p_value, method = "fdr"),
      
      # keep your normalization of overlap_genes, but make it safe
      overlap_genes = ifelse(
        is.na(overlap_genes) | overlap_genes == "",
        "",
        overlap_genes
      )
    ) %>%
    mutate(
      overlap_genes = ifelse(
        overlap_genes == "",
        "",
        overlap_genes %>%
          strsplit(",") %>%
          map(~ sort(.) %>% paste(collapse = ",")) %>%
          unlist()
      )
    )
  
  # convenience vectors (as before)
  original_overlap_genes <- unique(unlist(strsplit(original_df$overlap_genes, ",")))
  original_overlap_genes <- original_overlap_genes[original_overlap_genes != ""]
  original_rows <- unique(original_df$Var1)
  original_cols <- unique(original_df$Var2)
  
  # ---------- (D) Significant ----------
  significant_df <- original_df %>%
    filter(fdr < fdr_threshold, gene_overlap_count >= overlap_threshold)
  
  significant_overlap_genes <- unique(unlist(strsplit(significant_df$overlap_genes, ",")))
  significant_overlap_genes <- significant_overlap_genes[significant_overlap_genes != ""]
  significant_rows <- unique(significant_df$Var1)
  significant_cols <- unique(significant_df$Var2)
  
  # ---------- (E) Unique significant ----------
  significant_uniq_df <- significant_df %>%
    group_by(overlap_genes, Var2) %>%
    nest() %>%
    mutate(data = map(data, ~ .x %>% arrange(fdr) %>% slice(1))) %>%
    unnest(cols = data) %>%
    ungroup() %>%
    select(
      Var1, Var2,
      Var1_n_genes, Var2_n_genes,          # NEW ✅
      p_value, chi2, odds_ratio, log2_odds_ratio,
      gene_overlap_count, expected_overlap, # NEW ✅
      overlap_genes, fdr
    )
  
  significant_uniq_overlap_genes <- unique(unlist(strsplit(significant_uniq_df$overlap_genes, ",")))
  significant_uniq_overlap_genes <- significant_uniq_overlap_genes[significant_uniq_overlap_genes != ""]
  significant_uniq_rows <- unique(significant_uniq_df$Var1)
  significant_uniq_cols <- unique(significant_uniq_df$Var2)
  
  # ---------- (F) Filtered matrices (for plotting compatibility) ----------
  original_list <- filter_matrices(data)
  
  # ---------- (G) Output ----------
  output_list <- list(
    original_data = list(
      list = original_list,
      df = original_df,
      rows = original_rows,
      cols = original_cols,
      overlap_genes = original_overlap_genes
    ),
    significant_data = list(
      df = significant_df,
      rows = significant_rows,
      cols = significant_cols,
      overlap_genes = significant_overlap_genes
    ),
    significant_uniq_data = list(
      df = significant_uniq_df,
      rows = significant_uniq_rows,
      cols = significant_uniq_cols,
      overlap_genes = significant_uniq_overlap_genes
    ),
    gene_list_sizes = gene_list_sizes
  )
  
  return(output_list)
}
