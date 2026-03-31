url <- "https://raw.githubusercontent.com/MaayanLab/Enrichr-Viz-Appyter/master/Enrichr-Processed-Library-Storage/Clustered_Scatterplots/GWAS_Catalog_2025.csv"
gwasCatalog2025 <- read.csv(url, stringsAsFactors = FALSE)


gwasCatalog2025 %>% head %>% 
  select(-c(x,y,cluster)) %>% 
  mutate(term = str_replace_all(term, " ", "_")) %>% 
  mutate(genes = str_replace_all(genes, " ", "|"))
  
  
gwasCatalog2025 %>%
  select(-c(x, y, cluster)) %>%
  mutate(
    term  = str_replace_all(term, " ", "_"),
    genes = str_replace_all(genes, " ", "|"),
    genes = str_split(genes, "\\|"),
    genes = lapply(genes, unique),  # opcjonalnie usuwa duplikaty
    n_genes = lengths(genes)        # liczy geny w każdej liście
  ) %>% 
  filter(n_genes >= 40) %>% 
  select(-n_genes) %>% 
  { set_names(.$genes, .$term) } -> gwasCatalog2025_list


run_full_overlap_analysis(
  gene_lists = c(flat_allGrSignatures_17.10.2025[c("global_GR_genes_globalDown5TissuesDerivedCells",
                                                   "global_GR_genes_globalUp5TissuesDerivedCells",
                                                   "minusGlobalDown5TissuesDerivedCells_BloodCellsDown",
                                                   "minusGlobalDown5TissuesDerivedCells_LungCellsDown",
                                                   "minusGlobalDown5TissuesDerivedCells_NeuralCellsDown",
                                                   "minusGlobalUp5TissuesDerivedCells_BloodCellsUp",
                                                   "minusGlobalUp5TissuesDerivedCells_LungCellsUp",
                                                   "minusGlobalUp5TissuesDerivedCells_NeuralCellsUp")]
                 , gwasCatalog2025_list),
  total_genes = hgnc_symbols_vector_v110,
  # total_genes = all_genes,
  rows_to_filter = names(gwasCatalog2025_list),
  cols_to_filter = c("global_GR_genes_globalDown5TissuesDerivedCells",
                     "global_GR_genes_globalUp5TissuesDerivedCells",
                     "minusGlobalDown5TissuesDerivedCells_BloodCellsDown",
                     "minusGlobalDown5TissuesDerivedCells_LungCellsDown",
                     "minusGlobalDown5TissuesDerivedCells_NeuralCellsDown",
                     "minusGlobalUp5TissuesDerivedCells_BloodCellsUp",
                     "minusGlobalUp5TissuesDerivedCells_LungCellsUp",
                     "minusGlobalUp5TissuesDerivedCells_NeuralCellsUp"),
  plot_title_or = "",
  verbose = F
) -> tmp


tmp$processed$original_data$df %>% 
  filter(gene_overlap_count > 2, 
         p_value < 0.05) %>% 
  filter(Var2 == "minusGlobalUp5TissuesDerivedCells_BloodCellsUp")



perform_chi2_tests <- function(
    datasets,
    total_genes,
    verbose = FALSE,
    fdr_mode = c("none", "row", "col"),
    epsilon = 1e-6,             # minimalna korekta Haldane–Anscombe
    memory_limit_gb = 100       # 💾 domyślny limit pamięci RAM (GB)
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
  # 🔹 Progress bar
  # ------------------------------------------------------------
  if (!verbose) {
    if (!requireNamespace("progress", quietly = TRUE)) {
      stop("Please install the 'progress' package: install.packages('progress')")
    }
    pb <- progress::progress_bar$new(
      format = "⏳ Running chi² tests [:bar] :percent (:current/:total) ETA: :eta",
      total = n_combinations, clear = FALSE, width = 70
    )
  }
  
  # ------------------------------------------------------------
  # 🔹 Funkcja pomocnicza do sprawdzania pamięci
  # ------------------------------------------------------------
  check_memory <- function(limit_gb) {
    used_gb <- sum(gc()[, 2]) / 1024  # w przybliżeniu MB -> GB
    if (used_gb > limit_gb) {
      stop(
        sprintf("❌ Memory limit exceeded: %.2f GB used (limit = %.1f GB).",
                used_gb, limit_gb)
      )
    }
  }
  
  # ------------------------------------------------------------
  # 🔹 Pętla po parach
  # ------------------------------------------------------------
  results <- lapply(seq_along(combinations), function(idx) {
    pair <- combinations[[idx]]
    if (!verbose) pb$tick()
    
    # 🧠 sprawdzenie pamięci
    check_memory(memory_limit_gb)
    
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
