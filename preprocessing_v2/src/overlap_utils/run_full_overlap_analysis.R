# ##############################################################################
# ---- Main function: Run full overlap analysis (Chi2 + log2(OR)) ----
# ##############################################################################

run_full_overlap_analysis <- function(
    gene_lists,
    total_genes = hgnc_symbols_vector_v110,
    plot_title = "Full overlap log2(OR) heatmap",
    epsilon = 1e-6,
    fdr_mode = "row",
    overlap_threshold = 2,
    fdr_threshold = 0.05,
    drawing_overlap_threshold = 3,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color_scale_range = NULL,
    palette = c("navy", "white", "firebrick3"),
    verbose = TRUE,
    return_all = TRUE,
    ...
) {
  #' Run full overlap analysis and visualize results
  #'
  #' This function performs chi-square tests for all pairwise overlaps
  #' between gene lists, computes odds ratios (and log2(OR)),
  #' processes the results into structured tables, and visualizes
  #' the log2(OR) matrix as a heatmap.
  #'
  #' @param gene_lists Named list of character vectors with gene symbols.
  #' @param total_genes Background vector of all genes considered in analysis.
  #' @param plot_title Title for the resulting heatmap.
  #' @param epsilon Small value added to each cell (Haldane–Anscombe correction).
  #' @param fdr_mode "row" (default) or "col"; FDR adjustment mode.
  #' @param overlap_threshold Minimum number of shared genes for significance.
  #' @param fdr_threshold FDR cutoff for significant results.
  #' @param drawing_overlap_threshold Minimum overlap count shown in the heatmap.
  #' @param cluster_rows, cluster_cols Logical; perform hierarchical clustering.
  #' @param color_scale_range Numeric range for color limits (log2(OR)).
  #' @param palette Color palette for the heatmap fill gradient.
  #' @param verbose Logical; print progress messages.
  #' @param return_all Logical; if TRUE, return all intermediate results.
  #' @param ... Additional arguments passed to `heatmap_overlap_log2OR_ggplot()`.
  #'
  #' @return A list containing:
  #'   - plot: ggplot2 heatmap object
  #'   - overlap: list of matrices and stats
  #'   - processed: processed overlap results (from processing_overlap_results)
  #'   - df: combined data frame of results (if return_all = TRUE)
  #'
  #' @examples
  #' run_full_overlap_analysis(list(
  #'   ListA = c("FKBP5", "NR3C1", "PER1"),
  #'   ListB = c("FKBP5", "TSC22D3", "SGK1"),
  #'   ListC = c("PER1", "SGK1", "FKBP4")
  #' ),
  #' palette = c("navy", "white", "firebrick3"),
  #' axis_text_angle = 90,
  #' text_size_axis = 12)
  
  # --- Basic validation ---
  if (!is.list(gene_lists))
    stop("❌ gene_lists must be a named list of gene vectors")
  if (is.null(names(gene_lists)))
    stop("❌ gene_lists must have names for labeling")
  
  if (verbose)
    message("▶ Starting overlap analysis for ", length(gene_lists), " gene lists...")
  
  # ============================================================
  # 1️⃣ Run chi-square tests (compute overlap, p-values, OR, log2OR)
  # ============================================================
  chi2_results <- perform_chi2_tests(
    datasets = gene_lists,
    total_genes = total_genes,
    fdr_mode = fdr_mode,
    epsilon = epsilon,
    verbose = verbose
  )
  
  if (verbose)
    message("▶ Processed overlap data: ",
            nrow(processed$original_data$df), " total comparisons")
  
  
  # ============================================================
  # 2️⃣ Process overlap results into data frames
  # ============================================================
  processed <- processing_overlap_results(
    data = chi2_results,
    genes_list = gene_lists,
    rows_to_filter = rownames(chi2_results$p_value_matrix),
    cols_to_filter = colnames(chi2_results$p_value_matrix),
    fdr_threshold = fdr_threshold,
    overlap_threshold = overlap_threshold
  )
  
  if (verbose)
    message("▶ Processed overlap data: ",
            nrow(processed$original_data$df), " total comparisons")
  
  # ============================================================
  # 3️⃣ Draw heatmap based on log2(OR)
  # ============================================================
  heatmap_plot <- heatmap_overlap_log2OR_ggplot(
    data_list = processed,
    data_type = "original_data",
    title = plot_title,
    palette = palette,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    color_scale_range = color_scale_range,
    overlap_threshold = drawing_overlap_threshold,
    ...
  )
  
  if (verbose)
    message("✅ Heatmap created successfully!")
  
  # ============================================================
  # 4️⃣ Return structured results
  # ============================================================
  result <- list(
    plot = heatmap_plot,
    overlap = chi2_results,
    processed = processed
  )
  
  if (return_all)
    result$df <- processed$original_data$df
  
  return(result)
}


run_full_overlap_analysis <- function(
    gene_lists,
    total_genes,
    # 🔹 tytuły
    plot_title_or = "",
    plot_title_chi2 = "",
    
    # 🔹 klastrowanie i wygląd
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    text_size_axis = 12,
    show_dendrograms = FALSE,
    
    # 🔹 parametry OR
    palette_or = c("#07243e", "#f5a3a3", "#8b0000"),
    color_scale_range_or = c(-5, 5),
    text_contrast_range_or = c(-3, 4.9),
    p_thresholds_or = c(0.01, 0.0001),
    color_rects_or = c("#97C426", "#2F4603"),
    
    # 🔹 parametry χ²
    palette_chi2 = c("white", "#f79d00", "#c62a00"),
    color_scale_range_chi2 = c(0, 12),
    text_contrast_range_chi2 = c(0, 4),
    p_thresholds_chi2 = c(0.01, 0.0001),
    color_rects_chi2 = c("#97C426", "#2F4603"),
    
    # 🔹 inne
    triangle_mode = "upper",
    row_labels_map = NULL,
    col_labels_map = NULL
) {
  # ============================================================
  # 📦 Pakiety
  # ============================================================
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  # ============================================================
  # 🧮 1️⃣ Przygotowanie danych
  # ============================================================
  list_names <- names(gene_lists)
  message(glue::glue("▶ Starting overlap analysis for {length(gene_lists)} gene lists..."))
  
  # tworzymy wszystkie możliwe pary porównań
  pairs <- expand.grid(list_names, list_names, stringsAsFactors = FALSE)
  colnames(pairs) <- c("ListA", "ListB")
  
  results <- list()
  
  for (i in seq_len(nrow(pairs))) {
    a <- pairs$ListA[i]
    b <- pairs$ListB[i]
    genes_a <- gene_lists[[a]]
    genes_b <- gene_lists[[b]]
    
    overlap_genes <- intersect(genes_a, genes_b)
    n_overlap <- length(overlap_genes)
    n_a <- length(genes_a)
    n_b <- length(genes_b)
    n_total <- length(total_genes)
    
    # contingency table
    table_2x2 <- matrix(c(
      n_overlap,
      n_a - n_overlap,
      n_b - n_overlap,
      n_total - (n_a + n_b - n_overlap)
    ), nrow = 2)
    
    suppressWarnings({
      chi2 <- suppressWarnings(chisq.test(table_2x2, correct = FALSE))
    })
    
    p_val <- chi2$p.value
    or <- (n_overlap / (n_a - n_overlap)) / (n_b / (n_total - n_b))
    log2or <- ifelse(or > 0, log2(or), -20)
    chi2_value <- as.numeric(chi2$statistic)
    
    results[[i]] <- data.frame(
      ListA = a,
      ListB = b,
      overlap = n_overlap,
      p_value = p_val,
      OR = or,
      log2OR = log2or,
      chi2_value = chi2_value,
      stringsAsFactors = FALSE
    )
    
    message(glue::glue("Chi²: {a} vs {b} → overlap: {n_overlap} p = {format(p_val, digits=3)} OR = {format(or, digits=3)} log2(OR) = {format(log2or, digits=3)}"))
  }
  
  df <- bind_rows(results)
  message(glue::glue("▶ Processed overlap data: {nrow(df)} total comparisons"))
  
  # ============================================================
  # 🧱 2️⃣ Konwersja do macierzy
  # ============================================================
  make_matrix <- function(value_col) {
    df %>%
      select(ListA, ListB, {{ value_col }}) %>%
      pivot_wider(names_from = ListB, values_from = {{ value_col }}) %>%
      column_to_rownames("ListA") %>%
      as.matrix()
  }
  
  data_matrices <- list(
    log2_odds_ratio_matrix = make_matrix(log2OR),
    number_overlap_matrix = make_matrix(overlap),
    p_value_matrix = make_matrix(p_value),
    chi2_value_matrix = make_matrix(chi2_value)
  )
  
  # ============================================================
  # 🧩 3️⃣ Tworzenie listy wynikowej
  # ============================================================
  processed <- list(
    original_data = list(
      list = data_matrices,
      cols = colnames(data_matrices$log2_odds_ratio_matrix),
      rows = rownames(data_matrices$log2_odds_ratio_matrix)
    ),
    gene_list_sizes = sapply(gene_lists, length)
  )
  
  # ============================================================
  # 🎨 4️⃣ Generowanie wykresów
  # ============================================================
  plot_or <- heatmap_overlap_log2OR_ggplot(
    data_list = processed,
    data_type = "original_data",
    color_scale_range = color_scale_range_or,
    text_contrast_range = text_contrast_range_or,
    palette = palette_or,
    triangle_mode = triangle_mode,
    title = plot_title_or,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    text_size_axis = text_size_axis,
    show_dendrograms = show_dendrograms,
    row_labels_map = row_labels_map,
    col_labels_map = col_labels_map,
    p_thresholds = p_thresholds_or,
    color_rects = color_rects_chi2
  )
  
  plot_chi2 <- heatmap_overlap_log2CHI2_ggplot(
    data_list = processed,
    data_type = "original_data",
    color_scale_range = color_scale_range_chi2,
    text_contrast_range = text_contrast_range_chi2,
    palette = palette_chi2,
    triangle_mode = triangle_mode,
    title = plot_title_chi2,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    text_size_axis = text_size_axis,
    show_dendrograms = show_dendrograms,
    p_thresholds = p_thresholds_chi2,
    color_rects = color_rects_chi2,
    row_labels_map = row_labels_map,
    col_labels_map = col_labels_map,
  )
  
  # ============================================================
  # 📦 5️⃣ Zwracanie wyników
  # ============================================================
  return(list(
    processed = processed,
    plot_or = plot_or,
    plot_chi2 = plot_chi2
  ))
}


run_full_overlap_analysis <- function(
    gene_lists,
    total_genes,
    # 🔹 tytuły
    plot_title_or = "",
    plot_title_chi2 = "",
    
    # 🔹 klastrowanie i wygląd
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    text_size_axis = 12,
    show_dendrograms = FALSE,
    
    # 🔹 parametry OR
    palette_or = c("#07243e", "#f5a3a3", "#8b0000"),
    color_scale_range_or = c(-5, 5),
    text_contrast_range_or = c(-3, 4.9),
    p_thresholds_or = c(0.01, 0.0001),
    color_rects_or = c("#97C426", "#2F4603"),
    
    # 🔹 parametry χ²
    palette_chi2 = c("white", "#f79d00", "#c62a00"),
    color_scale_range_chi2 = c(0, 12),
    text_contrast_range_chi2 = c(0, 4),
    p_thresholds_chi2 = c(0.01, 0.0001),
    color_rects_chi2 = c("#97C426", "#2F4603"),
    
    # 🔹 inne
    fdr_threshold = 0.05,
    overlap_threshold = 2,
    triangle_mode = "upper",
    row_labels_map = NULL,
    col_labels_map = NULL
) {
  # ============================================================
  # 📦 Pakiety
  # ============================================================
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(glue)
  library(purrr)
  
  # ============================================================
  # 🧮 1️⃣ Obliczenia parowe
  # ============================================================
  list_names <- names(gene_lists)
  message(glue("▶ Starting overlap analysis for {length(gene_lists)} gene lists..."))
  
  pairs <- expand.grid(ListA = list_names, ListB = list_names, stringsAsFactors = FALSE)
  
  results <- purrr::pmap_dfr(pairs, function(ListA, ListB) {
    genes_a <- gene_lists[[ListA]]
    genes_b <- gene_lists[[ListB]]
    
    overlap_genes <- intersect(genes_a, genes_b)
    n_overlap <- length(overlap_genes)
    n_a <- length(genes_a)
    n_b <- length(genes_b)
    n_total <- length(total_genes)
    
    table_2x2 <- matrix(c(
      n_overlap,
      n_a - n_overlap,
      n_b - n_overlap,
      n_total - (n_a + n_b - n_overlap)
    ), nrow = 2)
    
    suppressWarnings({
      chi2 <- suppressWarnings(chisq.test(table_2x2, correct = FALSE))
    })
    
    p_val <- chi2$p.value
    or <- (n_overlap / (n_a - n_overlap + 1e-6)) / (n_b / (n_total - n_b + 1e-6))
    log2or <- ifelse(or > 0, log2(or), -20)
    chi2_value <- as.numeric(chi2$statistic)
    
    data.frame(
      ListA = ListA,
      ListB = ListB,
      overlap = n_overlap,
      p_value = p_val,
      OR = or,
      log2OR = log2or,
      chi2_value = chi2_value,
      overlap_genes = paste(sort(overlap_genes), collapse = ","),
      stringsAsFactors = FALSE
    )
  })
  
  message(glue("▶ Processed {nrow(results)} total comparisons"))
  
  # ============================================================
  # 🧩 2️⃣ Tworzenie macierzy
  # ============================================================
  make_matrix <- function(value_col) {
    results %>%
      select(ListA, ListB, {{ value_col }}) %>%
      pivot_wider(names_from = ListB, values_from = {{ value_col }}) %>%
      column_to_rownames("ListA") %>%
      as.matrix()
  }
  
  data_matrices <- list(
    log2_odds_ratio_matrix = make_matrix(log2OR),
    number_overlap_matrix  = make_matrix(overlap),
    p_value_matrix         = make_matrix(p_value),
    chi2_value_matrix      = make_matrix(chi2_value),
    overlap_genes_matrix   = make_matrix(overlap_genes)
  )
  
  # ============================================================
  # 🧮 3️⃣ Obliczenie FDR i filtrowanie
  # ============================================================
  results <- results %>%
    mutate(
      fdr = p.adjust(p_value, method = "fdr"),
      overlap_genes_char = ifelse(overlap_genes == "", NA, gsub(",", "|", overlap_genes))
    )
  
  # original_data
  original_df <- results
  original_overlap_genes <- unique(unlist(strsplit(original_df$overlap_genes, ",")))
  
  # significant_data
  significant_df <- original_df %>%
    filter(fdr < fdr_threshold, overlap >= overlap_threshold)
  significant_overlap_genes <- unique(unlist(strsplit(significant_df$overlap_genes, ",")))
  
  # significant unique (per Var2, top by fdr)
  significant_uniq_df <- significant_df %>%
    group_by(ListB, overlap_genes) %>%
    slice_min(fdr, n = 1) %>%
    ungroup()
  significant_uniq_overlap_genes <- unique(unlist(strsplit(significant_uniq_df$overlap_genes, ",")))
  
  # ============================================================
  # 🧱 4️⃣ Struktura wynikowa
  # ============================================================
  gene_list_sizes <- sapply(gene_lists, length)
  
  final_output <- list(
    original_data = list(
      df = original_df,
      rows = unique(original_df$ListA),
      cols = unique(original_df$ListB),
      overlap_genes = original_overlap_genes,
      list = data_matrices
    ),
    significant_data = list(
      df = significant_df,
      rows = unique(significant_df$ListA),
      cols = unique(significant_df$ListB),
      overlap_genes = significant_overlap_genes
    ),
    significant_uniq_data = list(
      df = significant_uniq_df,
      rows = unique(significant_uniq_df$ListA),
      cols = unique(significant_uniq_df$ListB),
      overlap_genes = significant_uniq_overlap_genes
    ),
    gene_list_sizes = gene_list_sizes
  )
  
  # ============================================================
  # 🎨 5️⃣ Generowanie wykresów
  # ============================================================
  plot_or <- heatmap_overlap_log2OR_ggplot(
    data_list = list(original_data = final_output$original_data),
    data_type = "original_data",
    color_scale_range = color_scale_range_or,
    text_contrast_range = text_contrast_range_or,
    palette = palette_or,
    triangle_mode = triangle_mode,
    title = plot_title_or,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    text_size_axis = text_size_axis,
    show_dendrograms = show_dendrograms,
    row_labels_map = row_labels_map,
    col_labels_map = col_labels_map,
    p_thresholds = p_thresholds_or,
    color_rects = color_rects_or
  )
  
  plot_chi2 <- heatmap_overlap_log2CHI2_ggplot(
    data_list = list(original_data = final_output$original_data),
    data_type = "original_data",
    color_scale_range = color_scale_range_chi2,
    text_contrast_range = text_contrast_range_chi2,
    palette = palette_chi2,
    triangle_mode = triangle_mode,
    title = plot_title_chi2,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    text_size_axis = text_size_axis,
    show_dendrograms = show_dendrograms,
    p_thresholds = p_thresholds_chi2,
    color_rects = color_rects_chi2,
    row_labels_map = row_labels_map,
    col_labels_map = col_labels_map
  )
  
  # ============================================================
  # 📦 6️⃣ Zwracanie wyników
  # ============================================================
  return(list(
    original_data = final_output$original_data,
    significant_data = final_output$significant_data,
    significant_uniq_data = final_output$significant_uniq_data,
    gene_list_sizes = final_output$gene_list_sizes,
    plot_or = plot_or,
    plot_chi2 = plot_chi2
  ))
}

run_full_overlap_analysis <- function(
    gene_lists,
    total_genes,
    # 🔹 tytuły
    plot_title_or = "",
    plot_title_chi2 = "",
    
    # 🔹 klastrowanie i wygląd
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    text_size_axis = 12,
    show_dendrograms = FALSE,
    
    # 🔹 parametry OR
    palette_or = c("#07243e", "#f5a3a3", "#8b0000"),
    color_scale_range_or = c(-5, 5),
    text_contrast_range_or = c(-3, 4.9),
    p_thresholds_or = c(0.01, 0.0001),
    color_rects_or = c("#97C426", "#2F4603"),
    
    # 🔹 parametry χ²
    palette_chi2 = c("white", "#f79d00", "#c62a00"),
    color_scale_range_chi2 = c(0, 12),
    text_contrast_range_chi2 = c(0, 4),
    p_thresholds_chi2 = c(0.05, 0.01),
    color_rects_chi2 = c("#97C426", "#2F4603"),
    
    # 🔹 inne
    triangle_mode = "upper",
    row_labels_map = NULL,
    col_labels_map = NULL,
    rows_to_filter = NULL,
    cols_to_filter = NULL,
    
    verbose = TRUE,
    fdr_mode = "row",
    epsilon = 0.000001,
    overlap_threshold = 2,
    fdr_threshold = 0.05,
    data_type = "original_data"
) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  if (!is.list(gene_lists))
    stop("❌ gene_lists must be a named list of gene vectors")
  if (is.null(names(gene_lists)))
    stop("❌ gene_lists must have names for labeling")
  
  if (verbose)
    message("▶ Starting overlap analysis for ", length(gene_lists), " gene lists...")
  
  # ============================================================
  # 1️⃣ Run chi-square tests
  # ============================================================
  chi2_results <- perform_chi2_tests(
    datasets = gene_lists,
    total_genes = total_genes,
    fdr_mode = fdr_mode,
    epsilon = epsilon,
    verbose = verbose
  )
  
  message("▶ chi2 resutls, done")
  
  
  # ============================================================
  # 2️⃣ Ustal wiersze i kolumny analizy
  # ============================================================
  if (is.null(rows_to_filter))
    rows_to_filter <- rownames(chi2_results$p_value_matrix)
  if (is.null(cols_to_filter))
    cols_to_filter <- colnames(chi2_results$p_value_matrix)
  
  if (verbose) {
    message("▶ Rows selected: ", paste(rows_to_filter, collapse = ", "))
    message("▶ Cols selected: ", paste(cols_to_filter, collapse = ", "))
  }
  
  # ============================================================
  # 3️⃣ Process overlap results into data frames
  # ============================================================
  processed <- processing_overlap_results(
    data = chi2_results,
    genes_list = gene_lists,
    rows_to_filter = rows_to_filter,
    cols_to_filter = cols_to_filter,
    fdr_threshold = fdr_threshold,
    overlap_threshold = overlap_threshold
  )
  

  message("▶ Processed overlap data: ",
          nrow(processed$original_data$df), " total comparisons")
  
  # ============================================================
  # 4️⃣ Wykresy
  # ============================================================
  plot_or <- heatmap_overlap_log2OR_ggplot(
    data_list = processed,
    data_type = data_type,
    color_scale_range = color_scale_range_or,
    text_contrast_range = text_contrast_range_or,
    palette = palette_or,
    triangle_mode = triangle_mode,
    title = plot_title_or,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    text_size_axis = text_size_axis,
    show_dendrograms = show_dendrograms,
    row_labels_map = row_labels_map,
    col_labels_map = col_labels_map,
    p_thresholds = p_thresholds_or,
    color_rects = color_rects_or
  )
  
  plot_chi2 <- heatmap_overlap_log2CHI2_ggplot(
    data_list = processed,
    data_type = data_type,
    color_scale_range = color_scale_range_chi2,
    text_contrast_range = text_contrast_range_chi2,
    palette = palette_chi2,
    triangle_mode = triangle_mode,
    title = plot_title_chi2,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    text_size_axis = text_size_axis,
    show_dendrograms = show_dendrograms,
    p_thresholds = p_thresholds_chi2,
    color_rects = color_rects_chi2,
    row_labels_map = row_labels_map,
    col_labels_map = col_labels_map
  )
  
  # ============================================================
  # 5️⃣ Wynik
  # ============================================================
  return(list(
    processed = processed,
    plot_or = plot_or,
    plot_chi2 = plot_chi2
  ))
}

