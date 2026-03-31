perform_chi2_tests <- function(
    datasets,
    total_genes,
    verbose = FALSE,
    fdr_mode = c("none", "row", "col"),
    epsilon = 1e-6   # minimalna korekta Haldane–Anscombe
) {
  # ------------------------------------------------------------
  # 📘 Ustawienia wstępne
  # ------------------------------------------------------------
  fdr_mode <- match.arg(fdr_mode)
  n <- length(datasets)
  dataset_lengths <- sapply(datasets, length)
  names_list <- names(datasets)
  
  if (is.null(names_list) || any(names_list == "")) {
    stop("All datasets must be named for matrix labeling.")
  }
  
  # ------------------------------------------------------------
  # 🔹 Inicjalizacja macierzy
  # ------------------------------------------------------------
  overlap_matrix         <- matrix(0, n, n, dimnames = list(names_list, names_list))
  p_value_matrix         <- matrix(0, n, n, dimnames = list(names_list, names_list))
  chi2_value_matrix      <- matrix(0, n, n, dimnames = list(names_list, names_list))
  odds_ratio_matrix      <- matrix(NA, n, n, dimnames = list(names_list, names_list))
  log2_odds_ratio_matrix <- matrix(NA, n, n, dimnames = list(names_list, names_list))
  fdr_matrix             <- matrix(NA, n, n, dimnames = list(names_list, names_list))
  overlap_genes_matrix   <- matrix("", n, n, dimnames = list(names_list, names_list))
  
  # ------------------------------------------------------------
  # 🔹 Kombinacje par do testu
  # ------------------------------------------------------------
  combinations <- combn(names_list, 2, simplify = FALSE)
  n_combinations <- length(combinations)
  
  # ------------------------------------------------------------
  # 🔹 Progress bar (tylko gdy verbose = FALSE)
  # ------------------------------------------------------------
  if (!verbose) {
    if (!requireNamespace("progress", quietly = TRUE)) {
      stop("Please install the 'progress' package: install.packages('progress')")
    }
    pb <- progress::progress_bar$new(
      format = "⏳ Running chi² tests [:bar] :current/:total (:percent) ETA: :eta",
      total = n_combinations, clear = FALSE, width = 70
    )
  }
  
  # ------------------------------------------------------------
  # 🔹 Pętla po parach
  # ------------------------------------------------------------
  results <- lapply(combinations, function(pair) {
    if (!verbose) pb$tick()
    
    i <- match(pair[1], names_list)
    j <- match(pair[2], names_list)
    
    overlapping_genes <- intersect(datasets[[i]], datasets[[j]])
    overlap_count <- length(overlapping_genes)
    external_genes_count <- length(total_genes) - length(unique(c(datasets[[i]], datasets[[j]])))
    
    a <- overlap_count
    b <- dataset_lengths[i] - overlap_count
    c <- dataset_lengths[j] - overlap_count
    d <- external_genes_count
    
    # korekta Haldane–Anscombe
    a2 <- a + epsilon; b2 <- b + epsilon; c2 <- c + epsilon; d2 <- d + epsilon
    
    matrix_chi2 <- matrix(c(a, b, c, d), nrow = 2)
    test_result <- suppressWarnings(chisq.test(matrix_chi2))
    odds_ratio <- (a2 * d2) / (b2 * c2)
    log2_or <- log2(odds_ratio)
    
    if (verbose) {
      cat(
        sprintf(
          "Chi²: %s vs %s → overlap: %d | p = %.3g | OR = %.3f | log2(OR) = %.3f\n",
          pair[1], pair[2], a, test_result$p.value, odds_ratio, log2_or
        )
      )
    }
    
    list(
      i = i,
      j = j,
      overlap_count = overlap_count,
      p_value = test_result$p.value,
      chi2_value = test_result$statistic,
      odds_ratio = odds_ratio,
      log2_or = log2_or,
      overlapping_genes = overlapping_genes
    )
  })
  
  # ------------------------------------------------------------
  # 🔹 Wypełnianie macierzy wynikami
  # ------------------------------------------------------------
  for (result in results) {
    i <- result$i
    j <- result$j
    overlap_genes_str <- if (length(result$overlapping_genes) > 0) {
      paste(sort(unique(result$overlapping_genes)), collapse = ",")
    } else {
      NA_character_
    }
    
    overlap_matrix[i, j] <- overlap_matrix[j, i] <- result$overlap_count
    p_value_matrix[i, j] <- p_value_matrix[j, i] <- result$p_value
    chi2_value_matrix[i, j] <- chi2_value_matrix[j, i] <- result$chi2_value
    odds_ratio_matrix[i, j] <- odds_ratio_matrix[j, i] <- result$odds_ratio
    log2_odds_ratio_matrix[i, j] <- log2_odds_ratio_matrix[j, i] <- result$log2_or
    overlap_genes_matrix[i, j] <- overlap_genes_matrix[j, i] <- overlap_genes_str
  }
  overlap_genes_matrix[is.na(overlap_genes_matrix)] <- ""
  
  # ------------------------------------------------------------
  # 🔹 Korekcja FDR
  # ------------------------------------------------------------
  if (fdr_mode != "none") {
    for (k in 1:n) {
      if (fdr_mode == "row") {
        p_row <- p_value_matrix[k, ]
        fdr_matrix[k, ] <- p.adjust(p_row, method = "fdr")
      } else if (fdr_mode == "col") {
        p_col <- p_value_matrix[, k]
        fdr_matrix[, k] <- p.adjust(p_col, method = "fdr")
      }
    }
  }
  
  # ------------------------------------------------------------
  # 🔹 Wynik
  # ------------------------------------------------------------
  return(list(
    number_overlap_matrix = overlap_matrix,
    p_value_matrix = p_value_matrix,
    fdr_matrix = fdr_matrix,
    chi2_value_matrix = chi2_value_matrix,
    odds_ratio_matrix = odds_ratio_matrix,
    log2_odds_ratio_matrix = log2_odds_ratio_matrix,
    overlap_genes_matrix = overlap_genes_matrix
  ))
}

# processing_overlap_results <- function(
#     data,
#     genes_list,
#     rows_to_filter,
#     cols_to_filter,
#     fdr_threshold = 0.05,
#     overlap_threshold = 2,
#     total_genes = NULL,     # <- potrzebne, by ewentualnie policzyć OR
#     epsilon = 1e-6          # <- pseudo-count do OR
# ) {
#   # ============================================================
#   # processing_overlap_results (enhanced)
#   # - now computes OR/log2OR if missing in `data`
#   # ============================================================
#   
#   library(dplyr)
#   library(tidyr)
#   library(reshape2)
#   library(purrr)
#   library(tibble)
#   
#   # ---------- (A) If OR/log2OR are missing, compute them ----------
#   if (is.null(data$odds_ratio_matrix) || is.null(data$log2_odds_ratio_matrix)) {
#     if (is.null(total_genes)) {
#       stop("`total_genes` is required to compute odds ratios when they are missing in `data`.")
#     }
#     # sizes
#     n_bg <- length(total_genes)
#     list_sizes <- sapply(genes_list, length)
#     
#     # assume square matrices with same dimnames:
#     overlap_mat <- data$number_overlap_matrix
#     rn <- rownames(overlap_mat); cn <- colnames(overlap_mat)
#     or_mat <- matrix(NA_real_, nrow = nrow(overlap_mat), ncol = ncol(overlap_mat),
#                      dimnames = list(rn, cn))
#     log2or_mat <- or_mat
#     
#     for (i in seq_along(rn)) {
#       for (j in seq_along(cn)) {
#         a <- overlap_mat[i, j]
#         ni <- list_sizes[rn[i]]
#         nj <- list_sizes[cn[j]]
#         # union = ni + nj - a
#         d <- n_bg - (ni + nj - a)
#         b <- ni - a
#         c <- nj - a
#         # odds ratio with epsilon
#         or <- ((a + epsilon) * (d + epsilon)) / ((b + epsilon) * (c + epsilon))
#         or_mat[i, j] <- or
#         log2or_mat[i, j] <- log2(or)
#       }
#     }
#     data$odds_ratio_matrix <- or_mat
#     data$log2_odds_ratio_matrix <- log2or_mat
#   }
#   
#   # ---------- (B) internal helpers ----------
#   safe_melt <- function(x, nm) {
#     if (is.null(x)) return(NULL)
#     x[rows_to_filter, cols_to_filter, drop = FALSE] %>%
#       melt() %>%
#       `colnames<-`(c("Var1", "Var2", nm))
#   }
#   
#   filter_matrices <- function(results) {
#     out <- list()
#     # keep only those that exist
#     keep <- intersect(
#       names(results),
#       c("p_value_matrix","chi2_value_matrix","number_overlap_matrix",
#         "overlap_genes_matrix","odds_ratio_matrix","log2_odds_ratio_matrix","fdr_matrix")
#     )
#     for (nm in keep) {
#       out[[nm]] <- results[[nm]][rows_to_filter, cols_to_filter, drop = FALSE]
#     }
#     out
#   }
#   
#   # ---------- (C) Build original_df ----------
#   pieces <- list(
#     safe_melt(data$p_value_matrix,         "p_value"),
#     safe_melt(data$chi2_value_matrix,      "chi2"),
#     safe_melt(data$number_overlap_matrix,  "gene_overlap_count"),
#     safe_melt(data$overlap_genes_matrix,   "overlap_genes"),
#     safe_melt(data$odds_ratio_matrix,      "odds_ratio"),
#     safe_melt(data$log2_odds_ratio_matrix, "log2_odds_ratio"),
#     safe_melt(data$fdr_matrix,             "fdr_value")
#   )
#   # drop NULLs and reduce
#   pieces <- pieces[!vapply(pieces, is.null, logical(1))]
#   original_df <- Reduce(function(x, y) dplyr::left_join(x, y, by = c("Var1","Var2")), pieces)
#   
#   # add computed FDR (row-wise p.adjust on the melted p-values)
#   original_df <- original_df %>%
#     mutate(
#       Var1 = as.character(Var1),
#       Var2 = as.character(Var2),
#       fdr = p.adjust(p_value, method = "fdr"),
#       overlap_genes = overlap_genes %>%
#         strsplit(",") %>%
#         map(~ sort(.) %>% paste(collapse = ",")) %>%
#         unlist()
#     )
#   
#   original_overlap_genes <- unique(unlist(strsplit(original_df$overlap_genes, ",")))
#   original_rows <- unique(original_df$Var1)
#   original_cols <- unique(original_df$Var2)
#   gene_list_sizes <- sapply(genes_list, length)
#   
#   # ---------- (D) Significant ----------
#   significant_df <- original_df %>%
#     filter(fdr < fdr_threshold, gene_overlap_count >= overlap_threshold)
#   
#   significant_overlap_genes <- unique(unlist(strsplit(significant_df$overlap_genes, ",")))
#   significant_rows <- unique(significant_df$Var1)
#   significant_cols <- unique(significant_df$Var2)
#   
#   # ---------- (E) Unique significant ----------
#   significant_uniq_df <- significant_df %>%
#     group_by(overlap_genes, Var2) %>%
#     nest() %>%
#     mutate(data = map(data, ~ .x %>% arrange(fdr) %>% slice(1))) %>%
#     unnest(cols = data) %>%
#     ungroup() %>%
#     select(Var1, Var2, p_value, chi2, odds_ratio, log2_odds_ratio,
#            gene_overlap_count, overlap_genes, fdr)
#   
#   significant_uniq_overlap_genes <- unique(unlist(strsplit(significant_uniq_df$overlap_genes, ",")))
#   significant_uniq_rows <- unique(significant_uniq_df$Var1)
#   significant_uniq_cols <- unique(significant_uniq_df$Var2)
#   
#   # ---------- (F) (Optional) provide filtered matrices for plotting compatibility ----------
#   original_list <- filter_matrices(data)
#   # (jeśli nie potrzebujesz macierzy w output — możesz usunąć to pole)
#   
#   # ---------- (G) Output ----------
#   output_list <- list(
#     original_data = list(
#       list = original_list,     # <- przydatne dla heatmap funkcji
#       df = original_df,
#       rows = original_rows,
#       cols = original_cols,
#       overlap_genes = original_overlap_genes
#     ),
#     significant_data = list(
#       df = significant_df,
#       rows = significant_rows,
#       cols = significant_cols,
#       overlap_genes = significant_overlap_genes
#     ),
#     significant_uniq_data = list(
#       df = significant_uniq_df,
#       rows = significant_uniq_rows,
#       cols = significant_uniq_cols,
#       overlap_genes = significant_uniq_overlap_genes
#     ),
#     gene_list_sizes = gene_list_sizes
#   )
#   
#   return(output_list)
# }
