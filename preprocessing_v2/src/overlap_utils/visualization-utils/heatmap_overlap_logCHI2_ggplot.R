heatmap_overlap_log2CHI2_ggplot <- function(
    data_list,
    data_type,
    p_thresholds = c(0.05, 0.01),
    color_rects = c("green", "red"),
    overlap_threshold = 3,
    col_only_n_genes = FALSE,
    col_mapping_vector = NULL,
    row_mapping_vector = NULL,
    gene_list_sizes = TRUE,
    col_significant = FALSE,
    color_scale_range = NULL,
    scale_color_limit = NULL,
    text_contrast_range = NULL,
    triangle_mode = c("full", "lower", "upper"),
    palette = c("navy", "white", "firebrick3"),
    text_color = "black",
    axis_text_angle = 45,
    axis_text_hjust = 0,
    text_size_tile = 3,
    text_size_axis = 10,
    text_size_legend = 10,
    title = NULL,
    title_size = 14,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_dendrograms = TRUE,
    cluster_method = "complete"
) {
  # ============================================================
  # 📦 Required packages
  # ============================================================
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggdendro)
  library(patchwork)
  
  triangle_mode <- match.arg(triangle_mode)
  
  # ============================================================
  # 1️⃣ Extract and validate input data
  # ============================================================
  data <- data_list[[data_type]]$list
  cols_to_filter <- data_list[[data_type]]$cols
  rows_to_filter <- data_list[[data_type]]$rows
  
  if (col_significant) {
    data <- lapply(data, function(mat) mat[rows_to_filter, cols_to_filter])
  }
  
  # --- Detect chi2 matrix name automatically ---
  if ("chi2_matrix" %in% names(data)) {
    chi2 <- as.matrix(data$chi2_matrix)
  } else if ("chi2_value_matrix" %in% names(data)) {
    chi2 <- as.matrix(data$chi2_value_matrix)
  } else {
    stop("❌ Missing chi2 matrix (expected: 'chi2_matrix' or 'chi2_value_matrix').")
  }
  
  # --- Check for other required matrices ---
  required <- c("log2_odds_ratio_matrix", "number_overlap_matrix", "p_value_matrix")
  missing <- required[!required %in% names(data)]
  if (length(missing) > 0) {
    stop(paste0(
      "❌ Missing required matrices: ", paste(missing, collapse = ", "),
      "\nExpected: log2_odds_ratio_matrix, number_overlap_matrix, p_value_matrix"
    ))
  }
  
  # --- Extract all matrices ---
  log2chi2 <- log2(chi2 + 1)
  log2or <- as.matrix(data$log2_odds_ratio_matrix)
  overlap <- as.matrix(data$number_overlap_matrix)
  pval <- as.matrix(data$p_value_matrix)
  
  # ============================================================
  # 2️⃣ Optional hierarchical clustering
  # ============================================================
  row_dend <- NULL
  col_dend <- NULL
  
  if (cluster_rows) {
    row_dend <- as.dendrogram(hclust(dist(log2chi2), method = cluster_method))
    rows_to_filter <- rownames(log2chi2)[order.dendrogram(row_dend)]
  }
  if (cluster_cols) {
    col_dend <- as.dendrogram(hclust(dist(t(log2chi2)), method = cluster_method))
    cols_to_filter <- colnames(log2chi2)[order.dendrogram(col_dend)]
  }
  
  # ============================================================
  # 3️⃣ Axis labels with gene list sizes
  # ============================================================
  if (gene_list_sizes) {
    if (is.null(col_mapping_vector)) {
      col_mapping_vector <- setNames(cols_to_filter, cols_to_filter)
    }
    col_mapping_vector <- setNames(
      if (col_only_n_genes) {
        as.character(data_list$gene_list_sizes[names(col_mapping_vector)])
      } else {
        paste0(col_mapping_vector, " (", data_list$gene_list_sizes[names(col_mapping_vector)], ")")
      },
      names(col_mapping_vector)
    )
    
    if (is.null(row_mapping_vector)) {
      row_mapping_vector <- setNames(rows_to_filter, rows_to_filter)
    }
    row_mapping_vector <- setNames(
      paste0(row_mapping_vector, " (", data_list$gene_list_sizes[names(row_mapping_vector)], ")"),
      names(row_mapping_vector)
    )
  }
  
  # ============================================================
  # 4️⃣ Prepare dataframe for ggplot
  # ============================================================
  df <- expand.grid(
    row = rownames(log2chi2),
    col = colnames(log2chi2),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      log2chi2 = as.vector(log2chi2),
      log2or = as.vector(log2or),
      overlap = as.vector(overlap),
      pval = as.vector(pval),
      row_label = row_mapping_vector[row],
      label = sprintf("%.2f (%d)", log2chi2, overlap)
    )
  
  df$row_factor <- factor(df$row, levels = rows_to_filter)
  df$col_factor <- factor(df$col, levels = cols_to_filter)
  
  # ============================================================
  # 5️⃣ Optional triangle filtering
  # ============================================================
  if (triangle_mode == "lower") {
    df <- df %>% filter(as.numeric(row_factor) > as.numeric(col_factor))
  } else if (triangle_mode == "upper") {
    df <- df %>% filter(as.numeric(row_factor) < as.numeric(col_factor))
  }
  
  # ============================================================
  # 6️⃣ Significance highlighting (p-value + log2OR)
  # ============================================================
  for (i in seq_along(p_thresholds)) {
    df[[paste0("sig_", i)]] <- (
      df$pval < p_thresholds[i] &
        df$overlap >= overlap_threshold &
        df$log2or > 0
    )
  }
  
  # ============================================================
  # 7️⃣ Color scaling and clipping
  # ============================================================
  if (is.null(color_scale_range)) {
    max_abs <- max(df$log2chi2, na.rm = TRUE)
    color_scale_range <- c(0, max_abs)
  }
  
  if (!is.null(scale_color_limit)) {
    df$log2chi2 <- pmax(pmin(df$log2chi2, scale_color_limit[2]), scale_color_limit[1])
    color_scale_range <- scale_color_limit
  }
  
  # ============================================================
  # 8️⃣ Dynamic text color (contrast range)
  # ============================================================
  if (!is.null(text_contrast_range)) {
    df$text_col <- ifelse(
      df$log2chi2 >= text_contrast_range[1] & df$log2chi2 <= text_contrast_range[2],
      "black",
      "white"
    )
  } else {
    df$text_col <- text_color
  }
  
  # ============================================================
  # 9️⃣ Build ggplot heatmap
  # ============================================================
  p_heatmap <- ggplot(df, aes(x = col_factor, y = row_factor)) +
    geom_tile(aes(fill = log2chi2), color = "white") +
    geom_text(aes(label = label, color = text_col), size = text_size_tile) +
    scale_color_identity()
  
  # --- Add significance borders
  for (i in seq_along(p_thresholds)) {
    sig_df <- df %>% filter(!!sym(paste0("sig_", i)))
    if (nrow(sig_df) > 0) {
      p_heatmap <- p_heatmap + geom_tile(
        data = sig_df,
        color = color_rects[i],
        fill = NA,
        size = 1.2
      )
    }
  }
  
  # ============================================================
  # 🔟 Final styling and theme
  # ============================================================
  p_heatmap <- p_heatmap +
    scale_fill_gradientn(colors = palette, limits = color_scale_range, name = "log2(χ² + 1)") +
    scale_x_discrete(position = "top", labels = col_mapping_vector[cols_to_filter]) +
    scale_y_discrete(position = "right", labels = row_mapping_vector[rows_to_filter]) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = axis_text_angle, hjust = axis_text_hjust, color = text_color, size = text_size_axis),
      axis.text.y = element_text(color = text_color, size = text_size_axis),
      legend.title = element_text(size = text_size_legend, color = text_color),
      legend.text = element_text(size = text_size_legend, color = text_color),
      panel.grid = element_blank(),
      axis.title = element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "vertical",
      plot.margin = margin(0, 0, 0, 0)
    )
  
  if (!is.null(title)) {
    p_heatmap <- p_heatmap + ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5, vjust = 1.2, size = title_size, color = text_color))
  }
  
  return(p_heatmap)
}

heatmap_overlap_log2CHI2_ggplot <- function(
    data_list,
    data_type,
    p_thresholds = c(0.05, 0.01),
    color_rects = c("green", "red"),
    overlap_threshold = 3,
    triangle_mode = c("full", "lower", "upper"),
    palette = c("white", "#f79d00", "#c62a00"),
    color_scale_range = NULL,
    text_contrast_range = NULL,
    text_color = "black",
    axis_text_angle = 45,
    axis_text_hjust = 0,
    text_size_tile = 3,
    text_size_axis = 10,
    text_size_legend = 10,
    title = NULL,
    title_size = 14,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_dendrograms = FALSE,
    cluster_method = "complete",
    # 👇 etykietowanie (tak jak w wersji OR)
    row_labels_map = NULL,
    col_labels_map = NULL
) {
  # ============================================================
  # 📦 Pakiety
  # ============================================================
  library(ggplot2)
  library(dplyr)
  library(scales)
  
  triangle_mode <- match.arg(triangle_mode)
  data <- data_list[[data_type]]$list
  
  # ============================================================
  # 🧮 Wczytanie danych
  # ============================================================
  if ("chi2_matrix" %in% names(data)) {
    chi2 <- as.matrix(data$chi2_matrix)
  } else if ("chi2_value_matrix" %in% names(data)) {
    chi2 <- as.matrix(data$chi2_value_matrix)
  } else {
    stop("❌ Missing chi2 matrix (expected: 'chi2_matrix' or 'chi2_value_matrix').")
  }
  
  required <- c("log2_odds_ratio_matrix", "number_overlap_matrix", "p_value_matrix")
  missing <- required[!required %in% names(data)]
  if (length(missing) > 0) {
    stop(paste0("❌ Missing required matrices: ", paste(missing, collapse = ", ")))
  }
  
  log2chi2 <- log2(chi2 + 1)
  log2or <- as.matrix(data$log2_odds_ratio_matrix)
  overlap <- as.matrix(data$number_overlap_matrix)
  pval <- as.matrix(data$p_value_matrix)
  
  # ============================================================
  # 🧬 Nazwy i rozmiary list genów
  # ============================================================
  original_rows <- rownames(log2chi2)
  original_cols <- colnames(log2chi2)
  gene_sizes <- data_list$gene_list_sizes
  
  # ============================================================
  # 🏷️ Mapowanie nazw (tak samo jak w OR)
  # ============================================================
  if (!is.null(row_labels_map)) {
    rownames(log2chi2) <- ifelse(original_rows %in% names(row_labels_map),
                                 row_labels_map[original_rows], original_rows)
  }
  if (!is.null(col_labels_map)) {
    colnames(log2chi2) <- ifelse(original_cols %in% names(col_labels_map),
                                 col_labels_map[original_cols], original_cols)
  }
  
  rows_to_filter <- rownames(log2chi2)
  cols_to_filter <- colnames(log2chi2)
  
  # ============================================================
  # 📊 Etykiety z liczbą genów
  # ============================================================
  safe_lookup <- function(x) gene_sizes[match(x, names(gene_sizes))]
  
  row_mapping_vector <- setNames(
    paste0(rows_to_filter, " (", safe_lookup(original_rows), ")"),
    rows_to_filter
  )
  col_mapping_vector <- setNames(
    paste0(cols_to_filter, " (", safe_lookup(original_cols), ")"),
    cols_to_filter
  )
  
  # ============================================================
  # 🔢 Dane do ggplot
  # ============================================================
  df <- expand.grid(
    row = rownames(log2chi2),
    col = colnames(log2chi2),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      log2chi2 = as.vector(log2chi2),
      log2or = as.vector(log2or),
      overlap = as.vector(overlap),
      pval = as.vector(pval),
      label = sprintf("%.2f (%d)", log2chi2, overlap)
    )
  
  df$row_factor <- factor(df$row, levels = rows_to_filter)
  df$col_factor <- factor(df$col, levels = cols_to_filter)
  
  # ============================================================
  # 🔺 Filtr trójkąta
  # ============================================================
  if (triangle_mode == "lower") {
    df <- df %>% filter(as.numeric(row_factor) > as.numeric(col_factor))
  } else if (triangle_mode == "upper") {
    df <- df %>% filter(as.numeric(row_factor) < as.numeric(col_factor))
  }
  
  # ============================================================
  # 📈 Oznaczanie istotnych pól
  # ============================================================
  for (i in seq_along(p_thresholds)) {
    df[[paste0("sig_", i)]] <- (
      df$pval < p_thresholds[i] &
        df$overlap >= overlap_threshold
    )
  }
  
  # ============================================================
  # 🎨 Skala kolorów
  # ============================================================
  if (is.null(color_scale_range)) {
    max_abs <- max(df$log2chi2, na.rm = TRUE)
    color_scale_range <- c(0, max_abs)
  }
  
  # przycinanie do zakresu
  df$log2chi2 <- pmax(pmin(df$log2chi2, color_scale_range[2]), color_scale_range[1])
  
  # ============================================================
  # 🔤 Kolor tekstu dynamicznie
  # ============================================================
  df$text_col <- if (!is.null(text_contrast_range)) {
    ifelse(df$log2chi2 >= text_contrast_range[1] & df$log2chi2 <= text_contrast_range[2],
           "black", "white")
  } else text_color
  
  # ============================================================
  # 🔥 Rysowanie heatmapy
  # ============================================================
  p <- ggplot(df, aes(x = col_factor, y = row_factor)) +
    geom_tile(aes(fill = log2chi2), color = "white") +
    geom_text(aes(label = label, color = text_col), size = text_size_tile) +
    scale_color_identity() +
    # 🔧 automatyczne przeskalowanie (rescale_mid tylko dla danych symetrycznych)
    scale_fill_gradientn(
      colors = palette,
      limits = color_scale_range,
      oob = scales::squish,
      rescaler = function(x, ...) {
        rng <- range(x, na.rm = TRUE)
        if (rng[1] < 0 && rng[2] > 0) {
          scales::rescale_mid(x, mid = 0)
        } else {
          scales::rescale(x, from = rng)
        }
      },
      na.value = tail(palette, 1),
      name = "log2(χ² + 1)"
    )
  
  # ============================================================
  # 🟥 Ramki dla istotnych pól
  # ============================================================
  for (i in seq_along(p_thresholds)) {
    sig_df <- df %>% filter(!!sym(paste0("sig_", i)))
    if (nrow(sig_df) > 0) {
      p <- p + geom_tile(
        data = sig_df,
        color = color_rects[i],
        fill = NA,
        size = 1.1
      )
    }
  }
  
  # ============================================================
  # 🧱 Styl i etykiety
  # ============================================================
  p <- p +
    scale_x_discrete(position = "top", labels = col_mapping_vector[cols_to_filter]) +
    scale_y_discrete(position = "right", labels = row_mapping_vector[rows_to_filter]) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = axis_text_angle, hjust = axis_text_hjust, color = text_color, size = text_size_axis),
      axis.text.y = element_text(color = text_color, size = text_size_axis),
      legend.title = element_text(size = text_size_legend, color = text_color),
      legend.text = element_text(size = text_size_legend, color = text_color),
      panel.grid = element_blank(),
      axis.title = element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "vertical"
    )
  
  if (!is.null(title)) {
    p <- p + ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5, vjust = 1.2, size = title_size, color = text_color))
  }
  
  return(p)
}


heatmap_overlap_log2CHI2_ggplot <- function(
    data_list,
    data_type,
    # 🔹 Nowe argumenty
    rows_to_filter = NULL,
    cols_to_filter = NULL,
    
    # 🔹 Istniejące parametry
    p_thresholds = c(0.05, 0.01),
    color_rects = c("green", "red"),
    overlap_threshold = 3,
    triangle_mode = c("full", "lower", "upper"),
    palette = c("white", "#f79d00", "#c62a00"),
    color_scale_range = NULL,
    text_contrast_range = NULL,
    text_color = "black",
    axis_text_angle = 45,
    axis_text_hjust = 0,
    text_size_tile = 3,
    text_size_axis = 10,
    text_size_legend = 10,
    title = NULL,
    title_size = 14,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_dendrograms = FALSE,
    cluster_method = "complete",
    row_labels_map = NULL,
    col_labels_map = NULL
) {
  # ============================================================
  # 📦 Pakiety
  # ============================================================
  library(ggplot2)
  library(dplyr)
  library(scales)
  
  triangle_mode <- match.arg(triangle_mode)
  data <- data_list[[data_type]]$list
  
  # ============================================================
  # 🧮 Pobranie macierzy
  # ============================================================
  chi2 <- if ("chi2_matrix" %in% names(data)) {
    as.matrix(data$chi2_matrix)
  } else if ("chi2_value_matrix" %in% names(data)) {
    as.matrix(data$chi2_value_matrix)
  } else stop("❌ Missing chi2 matrix.")
  
  log2or <- as.matrix(data$log2_odds_ratio_matrix)
  overlap <- as.matrix(data$number_overlap_matrix)
  pval <- as.matrix(data$p_value_matrix)
  
  # ============================================================
  # ✂️ Ręczne filtrowanie wierszy/kolumn
  # ============================================================
  if (!is.null(rows_to_filter)) {
    keep_rows <- intersect(rownames(chi2), rows_to_filter)
    chi2 <- chi2[keep_rows, , drop = FALSE]
    log2or <- log2or[keep_rows, , drop = FALSE]
    overlap <- overlap[keep_rows, , drop = FALSE]
    pval <- pval[keep_rows, , drop = FALSE]
  }
  if (!is.null(cols_to_filter)) {
    keep_cols <- intersect(colnames(chi2), cols_to_filter)
    chi2 <- chi2[, keep_cols, drop = FALSE]
    log2or <- log2or[, keep_cols, drop = FALSE]
    overlap <- overlap[, keep_cols, drop = FALSE]
    pval <- pval[, keep_cols, drop = FALSE]
  }
  
  message("✅ Przycięto dane do wymiarów: ", nrow(chi2), " × ", ncol(chi2))
  
  # ============================================================
  # 🧮 Przygotowanie danych
  # ============================================================
  log2chi2 <- log2(chi2 + 1)
  
  rows_to_filter <- rownames(log2chi2)
  cols_to_filter <- colnames(log2chi2)
  
  # ============================================================
  # 🔢 Dane do ggplot
  # ============================================================
  df <- expand.grid(
    row = rows_to_filter,
    col = cols_to_filter,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      log2chi2 = as.vector(log2chi2),
      log2or = as.vector(log2or),
      overlap = as.vector(overlap),
      pval = as.vector(pval),
      label = sprintf("%.2f (%d)", log2chi2, overlap)
    )
  
  df$row_factor <- factor(df$row, levels = rows_to_filter)
  df$col_factor <- factor(df$col, levels = cols_to_filter)
  
  # ============================================================
  # 🔺 Filtr trójkąta
  # ============================================================
  if (triangle_mode == "lower") {
    df <- df %>% filter(as.numeric(row_factor) > as.numeric(col_factor))
  } else if (triangle_mode == "upper") {
    df <- df %>% filter(as.numeric(row_factor) < as.numeric(col_factor))
  }
  
  # ============================================================
  # 📈 Oznaczanie istotnych pól
  # ============================================================
  for (i in seq_along(p_thresholds)) {
    df[[paste0("sig_", i)]] <- (
      df$pval < p_thresholds[i] &
        df$overlap >= overlap_threshold
    )
  }
  
  # ============================================================
  # 🎨 Skala kolorów
  # ============================================================
  if (is.null(color_scale_range)) {
    max_abs <- max(df$log2chi2, na.rm = TRUE)
    color_scale_range <- c(0, max_abs)
  }
  df$log2chi2 <- pmax(pmin(df$log2chi2, color_scale_range[2]), color_scale_range[1])
  
  # ============================================================
  # 🔤 Kolor tekstu
  # ============================================================
  df$text_col <- if (!is.null(text_contrast_range)) {
    ifelse(df$log2chi2 >= text_contrast_range[1] & df$log2chi2 <= text_contrast_range[2],
           "black", "white")
  } else text_color
  
  # ============================================================
  # 🔥 Rysowanie heatmapy
  # ============================================================
  p <- ggplot(df, aes(x = col_factor, y = row_factor)) +
    geom_tile(aes(fill = log2chi2), color = "white") +
    geom_text(aes(label = label, color = text_col), size = text_size_tile) +
    scale_color_identity() +
    scale_fill_gradientn(
      colors = palette,
      limits = color_scale_range,
      oob = scales::squish,
      name = "log2(χ² + 1)"
    )
  
  # ============================================================
  # 🟥 Ramki dla istotnych pól
  # ============================================================
  for (i in seq_along(p_thresholds)) {
    sig_df <- df %>% filter(!!sym(paste0("sig_", i)))
    if (nrow(sig_df) > 0) {
      p <- p + geom_tile(
        data = sig_df,
        color = color_rects[i],
        fill = NA,
        size = 1.1
      )
    }
  }
  
  # ============================================================
  # 🧱 Styl
  # ============================================================
  p <- p +
    theme_minimal() +
    scale_x_discrete(position = "top") +
    scale_y_discrete(position = "right") +
    theme(
      axis.text.x = element_text(angle = axis_text_angle, hjust = axis_text_hjust, size = text_size_axis),
      axis.text.y = element_text(size = text_size_axis),
      legend.title = element_text(size = text_size_legend),
      legend.text = element_text(size = text_size_legend),
      panel.grid = element_blank(),
      axis.title = element_blank(),
      legend.position = "bottom"
    )
  
  if (!is.null(title)) {
    p <- p + ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5, size = title_size))
  }
  
  return(p)
}

