processing_overlap_results <- function(
    data,
    genes_list,
    rows_to_filter,
    cols_to_filter,
    fdr_threshold = 0.05,
    overlap_threshold = 2,
    total_genes = NULL,
    epsilon = 1e-6
) {
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  
  # ---------- (A) Uzupełnij brakujące OR/log2OR ----------
  if (is.null(data$odds_ratio_matrix) || is.null(data$log2_odds_ratio_matrix)) {
    if (is.null(total_genes)) stop("`total_genes` wymagane do obliczenia OR.")
    n_bg <- length(total_genes)
    list_sizes <- sapply(genes_list, length)
    overlap_mat <- data$number_overlap_matrix
    rn <- rownames(overlap_mat); cn <- colnames(overlap_mat)
    or_mat <- matrix(NA_real_, nrow = nrow(overlap_mat), ncol = ncol(overlap_mat),
                     dimnames = list(rn, cn))
    log2or_mat <- or_mat
    for (i in seq_along(rn)) {
      for (j in seq_along(cn)) {
        a <- overlap_mat[i, j]
        ni <- list_sizes[rn[i]]
        nj <- list_sizes[cn[j]]
        d <- n_bg - (ni + nj - a)
        b <- ni - a
        c <- nj - a
        or <- ((a + epsilon) * (d + epsilon)) / ((b + epsilon) * (c + epsilon))
        or_mat[i, j] <- or
        log2or_mat[i, j] <- log2(or)
      }
    }
    data$odds_ratio_matrix <- or_mat
    data$log2_odds_ratio_matrix <- log2or_mat
  }
  
  # ---------- (B) bezpieczne melt ----------
  safe_melt <- function(x, nm) {
    if (is.null(x)) return(NULL)
    if (length(dim(x)) < 2) return(NULL)
    df <- as.data.frame(as.table(x[rows_to_filter, cols_to_filter, drop = FALSE]))
    if (ncol(df) < 3) return(NULL)
    colnames(df)[1:3] <- c("Var1", "Var2", nm)
    df
  }
  
  # ---------- (C) filtrowanie macierzy ----------
  filter_matrices <- function(results) {
    keep <- intersect(
      names(results),
      c("p_value_matrix","chi2_value_matrix","number_overlap_matrix",
        "overlap_genes_matrix","odds_ratio_matrix","log2_odds_ratio_matrix","fdr_matrix")
    )
    out <- list()
    for (nm in keep) {
      out[[nm]] <- results[[nm]][rows_to_filter, cols_to_filter, drop = FALSE]
    }
    out
  }
  
  # ---------- (D) budowanie df ----------
  pieces <- list(
    safe_melt(data$p_value_matrix, "p_value"),
    safe_melt(data$chi2_value_matrix, "chi2"),
    safe_melt(data$number_overlap_matrix, "gene_overlap_count"),
    safe_melt(data$overlap_genes_matrix, "overlap_genes"),
    safe_melt(data$odds_ratio_matrix, "odds_ratio"),
    safe_melt(data$log2_odds_ratio_matrix, "log2_odds_ratio"),
    safe_melt(data$fdr_matrix, "fdr_value")
  )
  pieces <- pieces[!vapply(pieces, is.null, logical(1))]
  pieces <- lapply(pieces, function(df) {
    if (!"Var1" %in% names(df) & "row" %in% names(df)) names(df)[names(df) == "row"] <- "Var1"
    if (!"Var2" %in% names(df) & "column" %in% names(df)) names(df)[names(df) == "column"] <- "Var2"
    df
  })
  pieces <- Filter(function(df) all(c("Var1", "Var2") %in% names(df)), pieces)
  if (length(pieces) == 0) stop("Brak ramek z Var1/Var2.")
  original_df <- Reduce(function(x, y) dplyr::left_join(x, y, by = c("Var1","Var2")), pieces)
  
  # ---------- (E) dodaj FDR ----------
  original_df <- original_df %>%
    mutate(
      Var1 = as.character(Var1),
      Var2 = as.character(Var2),
      fdr = p.adjust(p_value, method = "fdr"),
      overlap_genes = overlap_genes %>%
        strsplit(",") %>%
        map(~ sort(.) %>% paste(collapse = ",")) %>%
        unlist()
    )
  
  # ---------- (F) filtruj istotne ----------
  significant_df <- original_df %>%
    filter(fdr < fdr_threshold, gene_overlap_count >= overlap_threshold)
  
  # ---------- (G) jeśli brak wyników — zakończ ----------
  if (nrow(significant_df) == 0) {
    message("⚠️ Brak istotnych wyników po filtracji — pomijam.")
    return(NULL)
  }
  
  significant_overlap_genes <- unique(unlist(strsplit(significant_df$overlap_genes, ",")))
  significant_rows <- unique(significant_df$Var1)
  significant_cols <- unique(significant_df$Var2)
  
  significant_uniq_df <- significant_df %>%
    group_by(overlap_genes, Var2) %>%
    nest() %>%
    mutate(data = map(data, ~ .x %>% arrange(fdr) %>% slice(1))) %>%
    unnest(cols = data) %>%
    ungroup() %>%
    select(Var1, Var2, p_value, chi2, odds_ratio, log2_odds_ratio,
           gene_overlap_count, overlap_genes, fdr)
  
  significant_uniq_overlap_genes <- unique(unlist(strsplit(significant_uniq_df$overlap_genes, ",")))
  significant_uniq_rows <- unique(significant_uniq_df$Var1)
  significant_uniq_cols <- unique(significant_uniq_df$Var2)
  
  # ---------- (H) dane do heatmap ----------
  original_list <- filter_matrices(data)
  
  # ---------- (I) zwrot ----------
  list(
    original_data = list(
      list = original_list,
      df = original_df,
      rows = unique(original_df$Var1),
      cols = unique(original_df$Var2),
      overlap_genes = unique(unlist(strsplit(original_df$overlap_genes, ",")))
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
    gene_list_sizes = sapply(genes_list, length)
  )
}



processing_overlap_results <- function(
    data,
    genes_list,
    rows_to_filter,
    cols_to_filter,
    fdr_threshold = 0.05,
    p_value_threshold = 0.05,
    overlap_threshold = 2,
    use_fdr = TRUE,
    require_positive_log2or = TRUE,
    total_genes = NULL,
    epsilon = 1e-6
) {
  library(dplyr); library(tidyr); library(purrr); library(tibble)
  
  msg <- function(...) message(sprintf(...))
  
  # ---------- helpers ----------
  get_mat <- function(lst, nm) if (!is.null(lst[[nm]])) lst[[nm]] else NULL
  subset_mat <- function(mat, r, c) {
    if (is.null(mat)) return(NULL)
    rr <- intersect(rownames(mat), r)
    cc <- intersect(colnames(mat), c)
    mat[rr, cc, drop = FALSE]
  }
  melt_mat <- function(mat, value_name) {
    if (is.null(mat)) return(NULL)
    df <- as.data.frame(as.table(mat))
    colnames(df) <- c("Var1","Var2", value_name)
    df
  }
  
  # ---------- OR/log2OR if missing ----------
  if (is.null(data$odds_ratio_matrix) || is.null(data$log2_odds_ratio_matrix)) {
    if (is.null(total_genes)) stop("`total_genes` wymagane do OR/log2OR.")
    msg("⏳ Wyliczam OR/log2OR ...")
    n_bg <- length(total_genes)
    list_sizes <- sapply(genes_list, length)
    ov <- data$number_overlap_matrix
    or_mat <- matrix(NA_real_, nrow(ov), ncol(ov), dimnames = dimnames(ov))
    log2or_mat <- or_mat
    for (i in seq_len(nrow(ov))) for (j in seq_len(ncol(ov))) {
      a <- ov[i, j]
      ni <- list_sizes[rownames(ov)[i]]
      nj <- list_sizes[colnames(ov)[j]]
      d <- n_bg - (ni + nj - a)
      b <- ni - a
      c <- nj - a
      or <- ((a + epsilon) * (d + epsilon)) / ((b + epsilon) * (c + epsilon))
      or_mat[i, j] <- or
      log2or_mat[i, j] <- log2(or)
    }
    data$odds_ratio_matrix <- or_mat
    data$log2_odds_ratio_matrix <- log2or_mat
    msg("✅ OR/log2OR gotowe.")
  }
  
  # ---------- 1) TWARDY SUBSET do rows_to_filter × cols_to_filter ----------
  msg("▶ Twarde przycięcie do %d wierszy × %d kolumn.", length(rows_to_filter), length(cols_to_filter))
  pval_mat     <- subset_mat(get_mat(data, "p_value_matrix"),           rows_to_filter, cols_to_filter)
  fdr_mat      <- subset_mat(get_mat(data, "fdr_matrix"),               rows_to_filter, cols_to_filter)
  chi2_mat     <- subset_mat(get_mat(data, "chi2_value_matrix"),        rows_to_filter, cols_to_filter)
  n_overlap    <- subset_mat(get_mat(data, "number_overlap_matrix"),    rows_to_filter, cols_to_filter)
  ov_genes_mat <- subset_mat(get_mat(data, "overlap_genes_matrix"),     rows_to_filter, cols_to_filter)
  or_mat       <- subset_mat(get_mat(data, "odds_ratio_matrix"),        rows_to_filter, cols_to_filter)
  log2or_mat   <- subset_mat(get_mat(data, "log2_odds_ratio_matrix"),   rows_to_filter, cols_to_filter)
  
  # sanity: musimy mieć p-value, log2OR i count
  if (is.null(pval_mat)) stop("Brak p_value_matrix po przycięciu.")
  if (is.null(log2or_mat)) stop("Brak log2_odds_ratio_matrix po przycięciu.")
  if (is.null(n_overlap)) stop("Brak number_overlap_matrix po przycięciu.")
  
  # ---------- 2) MASKA ISTOTNOŚCI komórkowej ----------
  if (use_fdr) {
    if (is.null(fdr_mat)) {
      msg("⚠️ Brak fdr_matrix — obliczam FDR z p-value.")
      fdr_vec <- p.adjust(as.numeric(pval_mat), method = "fdr")
      fdr_mat <- matrix(fdr_vec, nrow = nrow(pval_mat), ncol = ncol(pval_mat),
                        dimnames = dimnames(pval_mat))
    }
    sig_p <- fdr_mat < fdr_threshold
  } else {
    sig_p <- pval_mat < p_value_threshold
  }
  sig_overlap <- n_overlap >= overlap_threshold
  sig_log2or  <- if (require_positive_log2or) (log2or_mat > 0) else matrix(TRUE, nrow(pval_mat), ncol(pval_mat))
  cell_mask   <- sig_p & sig_overlap & sig_log2or
  
  msg("▶ Komórek spełniających kryteria: %d / %d", sum(cell_mask, na.rm = TRUE), length(cell_mask))
  
  # ---------- 3) WYBÓR WIERSZY/KOLUMN z ≥1 znaczącą komórką ----------
  keep_rows <- rownames(pval_mat)[rowSums(cell_mask, na.rm = TRUE) > 0]
  keep_cols <- colnames(pval_mat)[colSums(cell_mask, na.rm = TRUE) > 0]
  msg("▶ Zachowane wiersze: %d, kolumny: %d (z ≥1 komórką istotną).", length(keep_rows), length(keep_cols))
  
  # jeśli nic nie przeszło – zwróć NULL
  if (length(keep_rows) == 0 || length(keep_cols) == 0) {
    msg("⚠️ Brak istotnych komórek — zwracam NULL.")
    return(NULL)
  }
  
  # ---------- 4) ZBUDUJ DF (najpierw FULL w obrębie twardego subsetu), potem PRUNE ----------
  build_df <- function(pval, chi2, n_ov, ov_genes, orv, l2or, fdrm = NULL) {
    parts <- list(
      melt_mat(pval, "p_value"),
      melt_mat(chi2, "chi2"),
      melt_mat(n_ov, "gene_overlap_count"),
      melt_mat(ov_genes, "overlap_genes"),
      melt_mat(orv, "odds_ratio"),
      melt_mat(l2or, "log2_odds_ratio"),
      melt_mat(fdrm, "fdr_value")
    )
    parts <- parts[!vapply(parts, is.null, logical(1))]
    df <- Reduce(function(x, y) dplyr::left_join(x, y, by = c("Var1","Var2")), parts)
    df <- df %>%
      mutate(
        fdr = if (!is.null(fdrm)) fdr_value else p.adjust(p_value, method = "fdr"),
        overlap_genes = overlap_genes %||% "",
        overlap_genes = ifelse(is.na(overlap_genes), "", overlap_genes),
        overlap_genes = map_chr(strsplit(overlap_genes, ","), ~ paste(sort(.x[.x != ""]), collapse=",")),
        sig_1 = if (use_fdr) fdr < fdr_threshold else p_value < p_value_threshold
      )
    df
  }
  
  `%||%` <- function(a, b) if (is.null(a)) b else a
  
  # full (po twardym subset) i pruned (po masce)
  original_df <- build_df(pval_mat, chi2_mat, n_overlap, ov_genes_mat, or_mat, log2or_mat, fdr_mat)
  pruned_idx  <- paste(rep(keep_rows, each = length(keep_cols)), keep_cols, sep=".")
  original_df$._key <- paste(original_df$Var1, original_df$Var2, sep=".")
  significant_df <- original_df %>% filter(._key %in% pruned_idx) %>% select(-._key)
  
  # unikalne: po (overlap_genes, Var2) wybierz najniższy FDR (w pruned space)
  significant_uniq_df <- significant_df %>%
    group_by(overlap_genes, Var2) %>%
    arrange(fdr, .by_group = TRUE) %>%
    slice(1L) %>%
    ungroup()
  
  # ---------- 5) ZBUDUJ LISTY MACIERZY: original (subset R×C), significant (pruned R×C), uniq (pruned uniq R×C) ----------
  original_list <- list(
    p_value_matrix           = pval_mat,
    chi2_value_matrix        = chi2_mat,
    number_overlap_matrix    = n_overlap,
    overlap_genes_matrix     = ov_genes_mat,
    odds_ratio_matrix        = or_mat,
    log2_odds_ratio_matrix   = log2or_mat,
    fdr_matrix               = if (use_fdr) (if (!is.null(fdr_mat)) fdr_mat else {
      # oblicz z pval_mat
      matrix(p.adjust(as.numeric(pval_mat), method = "fdr"),
             nrow = nrow(pval_mat), ncol = ncol(pval_mat),
             dimnames = dimnames(pval_mat))
    }) else NULL
  )
  
  # przycięte macierze do keep_rows/keep_cols
  prune_list <- function(lst, r, c) {
    lapply(lst, function(m) if (is.null(m)) NULL else m[r, c, drop = FALSE])
  }
  significant_list      <- prune_list(original_list, keep_rows, keep_cols)
  
  # dla uniq bierz wiersze/kolumny faktycznie obecne po uniq
  uniq_rows <- unique(significant_uniq_df$Var1)
  uniq_cols <- unique(significant_uniq_df$Var2)
  significant_uniq_list <- prune_list(original_list, uniq_rows, uniq_cols)
  
  msg("✅ Finalne wymiary macierzy:")
  msg("   • original:           %d × %d", nrow(pval_mat), ncol(pval_mat))
  msg("   • significant(pruned): %d × %d", length(keep_rows), length(keep_cols))
  msg("   • significant_uniq:    %d × %d", length(uniq_rows), length(uniq_cols))
  
  # ---------- 6) ZWROT ----------
  list(
    original_data = list(
      list = original_list,
      df   = original_df,
      rows = rownames(pval_mat),
      cols = colnames(pval_mat),
      overlap_genes = unique(unlist(strsplit(original_df$overlap_genes, ",")))
    ),
    significant_data = list(
      list = significant_list,
      df   = significant_df,
      rows = keep_rows,
      cols = keep_cols,
      overlap_genes = unique(unlist(strsplit(significant_df$overlap_genes, ","))),
      filter_type = paste0(
        if (use_fdr) sprintf("FDR<%g", fdr_threshold) else sprintf("p<%g", p_value_threshold),
        if (require_positive_log2or) " & log2OR>0" else "",
        sprintf(" & overlap>=%d", overlap_threshold)
      )
    ),
    significant_uniq_data = list(
      list = significant_uniq_list,
      df   = significant_uniq_df,
      rows = uniq_rows,
      cols = uniq_cols,
      overlap_genes = unique(unlist(strsplit(significant_uniq_df$overlap_genes, ",")))
    ),
    gene_list_sizes = sapply(genes_list, length)
  )
}
