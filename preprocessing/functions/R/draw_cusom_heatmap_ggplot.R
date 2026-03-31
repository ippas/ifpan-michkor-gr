draw_custom_heatmap_ggplot_v2 <- function(data_list,
                                          data_type,
                                          p_thresholds = c(0.05, 0.01),
                                          color_rects = c("green", "red"),
                                          overlap_threshold = 3,
                                          apply_filling = TRUE,
                                          color_filling = "green",
                                          size_filling = 1,
                                          alpha_filling = 1,
                                          col_only_n_genes = FALSE,
                                          col_mapping_vector = NULL,
                                          row_mapping_vector = NULL,
                                          gene_list_sizes = TRUE,
                                          col_significant = FALSE,
                                          color_scale_range = NULL,
                                          palette = c("white", "#f1bcbb", "#edacab", "#e68a89")) {
  
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  
  data <- data_list[[data_type]]$list
  cols_to_filter <- data_list[[data_type]]$cols
  rows_to_filter <- data_list[[data_type]]$rows
  
  if (col_significant) {
    data <- lapply(data, function(mat) mat[rows_to_filter, cols_to_filter])
  }
  
  chi2 <- data$chi2_value_matrix
  overlap <- data$number_overlap_matrix
  pval <- data$p_value_matrix
  
  # Mapping
  if (gene_list_sizes) {
    if (is.null(col_mapping_vector)) {
      col_mapping_vector <- setNames(colnames(chi2), colnames(chi2))
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
      row_mapping_vector <- setNames(rownames(chi2), rownames(chi2))
    }
    row_mapping_vector <- setNames(
      paste0(row_mapping_vector, " (", data_list$gene_list_sizes[names(row_mapping_vector)], ")"),
      names(row_mapping_vector)
    )
  }
  
  # Prepare df
  df <- expand.grid(row = rownames(chi2), col = colnames(chi2), stringsAsFactors = FALSE) %>%
    mutate(
      chi2 = as.vector(chi2),
      overlap = as.vector(overlap),
      pval = as.vector(pval),
      row_label = row_mapping_vector[row],
      col_label = col_mapping_vector[col],
      label = sprintf("%.1f (%d)", chi2, overlap),
      log_chi2 = log2(chi2 + 1)
    )
  
  for (i in seq_along(p_thresholds)) {
    df[[paste0("sig_", i)]] <- df$pval < p_thresholds[i] & df$overlap >= overlap_threshold
  }
  
  if (is.null(color_scale_range)) {
    color_scale_range <- range(df$log_chi2, na.rm = TRUE)
  }
  
  # Black border coordinates (entire heatmap frame)
  border_df <- data.frame(
    x = c(min(df$col_label), max(df$col_label)),
    y = c(min(df$row_label), max(df$row_label))
  )
  
  p <- ggplot(df, aes(x = col_label, y = row_label)) +
    geom_tile(aes(fill = log_chi2), color = "white") +
    geom_text(aes(label = label), size = 3)
  
  for (i in seq_along(p_thresholds)) {
    sig_df <- df %>% filter(!!sym(paste0("sig_", i)))
    if (nrow(sig_df) > 0) {
      p <- p + geom_tile(data = sig_df, color = color_rects[i], fill = NA, size = 1.2)
    }
  }
  
  if (apply_filling) {
    p <- p + geom_point(data = df %>% filter(sig_1), shape = 16, size = size_filling,
                        color = color_filling, alpha = alpha_filling)
  }
  
  # Black frame
  p <- p + geom_rect(data = NULL,
                     aes(xmin = 0.5, xmax = length(unique(df$col_label)) + 0.5,
                         ymin = 0.5, ymax = length(unique(df$row_label)) + 0.5),
                     inherit.aes = FALSE,
                     color = "black", fill = NA, size = 1)
  
  p <- p +
    scale_fill_gradientn(colors = palette,
                         limits = color_scale_range,
                         name = "log2(chi2 + 1)") +
    scale_x_discrete(position = "top") +
    scale_y_discrete(position = "right") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 0),
      panel.grid = element_blank(),
      axis.title = element_blank()
    )
  
  return(p)
}

draw_custom_heatmap_ggplot_v2 <- function(data_list,
                                          data_type,
                                          p_thresholds = c(0.05, 0.01),
                                          color_rects = c("green", "red"),
                                          overlap_threshold = 3,
                                          apply_filling = TRUE,
                                          color_filling = "green",
                                          size_filling = 1,
                                          alpha_filling = 1,
                                          col_only_n_genes = FALSE,
                                          col_mapping_vector = NULL,
                                          row_mapping_vector = NULL,
                                          gene_list_sizes = TRUE,
                                          col_significant = FALSE,
                                          color_scale_range = NULL,
                                          palette = c("white", "#f1bcbb", "#edacab", "#e68a89")) {
  
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  
  data <- data_list[[data_type]]$list
  cols_to_filter <- data_list[[data_type]]$cols
  rows_to_filter <- data_list[[data_type]]$rows
  
  if (col_significant) {
    data <- lapply(data, function(mat) mat[rows_to_filter, cols_to_filter])
  }
  
  chi2 <- data$chi2_value_matrix
  overlap <- data$number_overlap_matrix
  pval <- data$p_value_matrix
  
  # Mapping
  if (gene_list_sizes) {
    if (is.null(col_mapping_vector)) {
      col_mapping_vector <- setNames(colnames(chi2), colnames(chi2))
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
      row_mapping_vector <- setNames(rownames(chi2), rownames(chi2))
    }
    row_mapping_vector <- setNames(
      paste0(row_mapping_vector, " (", data_list$gene_list_sizes[names(row_mapping_vector)], ")"),
      names(row_mapping_vector)
    )
  }
  
  # Prepare df
  df <- expand.grid(row = rownames(chi2), col = colnames(chi2), stringsAsFactors = FALSE) %>%
    mutate(
      chi2 = as.vector(chi2),
      overlap = as.vector(overlap),
      pval = as.vector(pval),
      row_label = row_mapping_vector[row],
      col_label = col_mapping_vector[col],
      label = sprintf("%.1f (%d)", chi2, overlap),
      log_chi2 = log2(chi2 + 1)
    )
  
  for (i in seq_along(p_thresholds)) {
    df[[paste0("sig_", i)]] <- df$pval < p_thresholds[i] & df$overlap >= overlap_threshold
  }
  
  if (is.null(color_scale_range)) {
    color_scale_range <- range(df$log_chi2, na.rm = TRUE)
  }
  
  # Start ggplot
  p <- ggplot(df, aes(x = col_label, y = row_label)) +
    geom_tile(aes(fill = log_chi2), color = "white") +
    geom_text(aes(label = label), size = 3)
  
  for (i in seq_along(p_thresholds)) {
    sig_df <- df %>% filter(!!sym(paste0("sig_", i)))
    if (nrow(sig_df) > 0) {
      p <- p + geom_tile(data = sig_df, color = color_rects[i], fill = NA, size = 1.2)
    }
  }
  
  if (apply_filling) {
    p <- p + geom_point(data = df %>% filter(sig_1), shape = 16, size = size_filling,
                        color = color_filling, alpha = alpha_filling)
  }
  
  # Black frame
  p <- p + geom_rect(data = NULL,
                     aes(xmin = 0.5, xmax = length(unique(df$col_label)) + 0.5,
                         ymin = 0.5, ymax = length(unique(df$row_label)) + 0.5),
                     inherit.aes = FALSE,
                     color = "black", fill = NA, size = 1)
  
  # Add scale and formatting
  p <- p +
    scale_fill_gradientn(colors = palette,
                         limits = color_scale_range,
                         name = "log2(chi2 + 1)") +
    scale_x_discrete(position = "top") +
    scale_y_discrete(position = "right") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 0),
      panel.grid = element_blank(),
      axis.title = element_blank(),
      legend.position = "bottom",          # <-- legenda POD spodem
      legend.direction = "horizontal",     # <-- poziomo
      legend.box = "vertical"              # <-- dla lepszego rozmieszczenia
    )
  
  return(p)
}

draw_custom_heatmap_ggplot_v3 <- function(data_list,
                                          data_type,
                                          p_thresholds = c(0.05, 0.01),
                                          color_rects = c("green", "red"),
                                          overlap_threshold = 3,
                                          apply_filling = TRUE,
                                          color_filling = "green",
                                          size_filling = 1,
                                          alpha_filling = 1,
                                          col_only_n_genes = FALSE,
                                          col_mapping_vector = NULL,
                                          row_mapping_vector = NULL,
                                          gene_list_sizes = TRUE,
                                          col_significant = FALSE,
                                          color_scale_range = NULL,
                                          palette = c("white", "#f1bcbb", "#edacab", "#e68a89"),
                                          # Nowe parametry:
                                          text_color = "black",
                                          axis_text_angle = 45,
                                          axis_text_hjust = 0,
                                          text_size_tile = 3,
                                          text_size_axis = 10,
                                          text_size_legend = 10) {
  
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  
  data <- data_list[[data_type]]$list
  cols_to_filter <- data_list[[data_type]]$cols
  rows_to_filter <- data_list[[data_type]]$rows
  
  if (col_significant) {
    data <- lapply(data, function(mat) mat[rows_to_filter, cols_to_filter])
  }
  
  chi2 <- data$chi2_value_matrix
  overlap <- data$number_overlap_matrix
  pval <- data$p_value_matrix
  
  # Mapping
  if (gene_list_sizes) {
    if (is.null(col_mapping_vector)) {
      col_mapping_vector <- setNames(colnames(chi2), colnames(chi2))
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
      row_mapping_vector <- setNames(rownames(chi2), rownames(chi2))
    }
    row_mapping_vector <- setNames(
      paste0(row_mapping_vector, " (", data_list$gene_list_sizes[names(row_mapping_vector)], ")"),
      names(row_mapping_vector)
    )
  }
  
  # Prepare df
  df <- expand.grid(row = rownames(chi2), col = colnames(chi2), stringsAsFactors = FALSE) %>%
    mutate(
      chi2 = as.vector(chi2),
      overlap = as.vector(overlap),
      pval = as.vector(pval),
      row_label = row_mapping_vector[row],
      col_label = col_mapping_vector[col],
      label = sprintf("%.1f (%d)", chi2, overlap),
      log_chi2 = log2(chi2 + 1)
    )
  
  for (i in seq_along(p_thresholds)) {
    df[[paste0("sig_", i)]] <- df$pval < p_thresholds[i] & df$overlap >= overlap_threshold
  }
  
  if (is.null(color_scale_range)) {
    color_scale_range <- range(df$log_chi2, na.rm = TRUE)
  }
  
  # Start ggplot
  p <- ggplot(df, aes(x = col_label, y = row_label)) +
    geom_tile(aes(fill = log_chi2), color = "white") +
    geom_text(aes(label = label), size = text_size_tile, color = text_color)
  
  for (i in seq_along(p_thresholds)) {
    sig_df <- df %>% filter(!!sym(paste0("sig_", i)))
    if (nrow(sig_df) > 0) {
      p <- p + geom_tile(data = sig_df, color = color_rects[i], fill = NA, size = 1.2)
    }
  }
  
  if (apply_filling) {
    p <- p + geom_point(data = df %>% filter(sig_1), shape = 16, size = size_filling,
                        color = color_filling, alpha = alpha_filling)
  }
  
  # Black frame
  p <- p + geom_rect(data = NULL,
                     aes(xmin = 0.5, xmax = length(unique(df$col_label)) + 0.5,
                         ymin = 0.5, ymax = length(unique(df$row_label)) + 0.5),
                     inherit.aes = FALSE,
                     color = "black", fill = NA, size = 1)
  
  # Final formatting
  p <- p +
    scale_fill_gradientn(colors = palette,
                         limits = color_scale_range,
                         name = "log2(chi2 + 1)") +
    scale_x_discrete(position = "top") +
    scale_y_discrete(position = "right") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = axis_text_angle, hjust = axis_text_hjust,
                                 color = text_color, size = text_size_axis),
      axis.text.y = element_text(color = text_color, size = text_size_axis),
      legend.title = element_text(size = text_size_legend, color = text_color),
      legend.text = element_text(size = text_size_legend, color = text_color),
      panel.grid = element_blank(),
      axis.title = element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "vertical"
    )
  
  return(p)
}

draw_custom_heatmap_ggplot_v3 <- function(data_list,
                                          data_type,
                                          p_thresholds = c(0.05, 0.01),
                                          color_rects = c("green", "red"),
                                          overlap_threshold = 3,
                                          apply_filling = TRUE,
                                          color_filling = "green",
                                          size_filling = 1,
                                          alpha_filling = 1,
                                          col_only_n_genes = FALSE,
                                          col_mapping_vector = NULL,
                                          row_mapping_vector = NULL,
                                          gene_list_sizes = TRUE,
                                          col_significant = FALSE,
                                          color_scale_range = NULL,
                                          palette = c("white", "#f1bcbb", "#edacab", "#e68a89"),
                                          # Nowe parametry:
                                          text_color = "black",
                                          axis_text_angle = 45,
                                          axis_text_hjust = 0,
                                          text_size_tile = 3,
                                          text_size_axis = 10,
                                          text_size_legend = 10) {
  
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  
  data <- data_list[[data_type]]$list
  cols_to_filter <- data_list[[data_type]]$cols
  rows_to_filter <- data_list[[data_type]]$rows
  
  if (col_significant) {
    data <- lapply(data, function(mat) mat[rows_to_filter, cols_to_filter])
  }
  
  chi2 <- data$chi2_value_matrix
  overlap <- data$number_overlap_matrix
  pval <- data$p_value_matrix
  
  # Mapping
  if (gene_list_sizes) {
    if (is.null(col_mapping_vector)) {
      col_mapping_vector <- setNames(colnames(chi2), colnames(chi2))
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
      row_mapping_vector <- setNames(rownames(chi2), rownames(chi2))
    }
    row_mapping_vector <- setNames(
      paste0(row_mapping_vector, " (", data_list$gene_list_sizes[names(row_mapping_vector)], ")"),
      names(row_mapping_vector)
    )
  }
  
  # Prepare df
  df <- expand.grid(row = rownames(chi2), col = colnames(chi2), stringsAsFactors = FALSE) %>%
    mutate(
      chi2 = as.vector(chi2),
      overlap = as.vector(overlap),
      pval = as.vector(pval),
      row_label = row_mapping_vector[row],
      col_label = col_mapping_vector[col],
      label = sprintf("%.1f (%d)", chi2, overlap),
      log_chi2 = log2(chi2 + 1)
    )
  
  for (i in seq_along(p_thresholds)) {
    df[[paste0("sig_", i)]] <- df$pval < p_thresholds[i] & df$overlap >= overlap_threshold
  }
  
  if (is.null(color_scale_range)) {
    color_scale_range <- range(df$log_chi2, na.rm = TRUE)
  }
  
  # Start ggplot
  p <- ggplot(df, aes(x = col_label, y = row_label)) +
    geom_tile(aes(fill = log_chi2), color = "white") +
    geom_text(aes(label = label), size = text_size_tile, color = text_color)
  
  for (i in seq_along(p_thresholds)) {
    sig_df <- df %>% filter(!!sym(paste0("sig_", i)))
    if (nrow(sig_df) > 0) {
      p <- p + geom_tile(data = sig_df, color = color_rects[i], fill = NA, size = 1.2)
    }
  }
  
  if (apply_filling) {
    p <- p + geom_point(data = df %>% filter(sig_1), shape = 16, size = size_filling,
                        color = color_filling, alpha = alpha_filling)
  }
  
  # Black frame
  p <- p + geom_rect(data = NULL,
                     aes(xmin = 0.5, xmax = length(unique(df$col_label)) + 0.5,
                         ymin = 0.5, ymax = length(unique(df$row_label)) + 0.5),
                     inherit.aes = FALSE,
                     color = "black", fill = NA, size = 1)
  
  # Final formatting
  p <- p +
    scale_fill_gradientn(colors = palette,
                         limits = color_scale_range,
                         name = "log2(chi2 + 1)") +
    scale_x_discrete(position = "top") +
    scale_y_discrete(position = "right") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = axis_text_angle, hjust = axis_text_hjust,
                                 color = text_color, size = text_size_axis),
      axis.text.y = element_text(color = text_color, size = text_size_axis),
      legend.title = element_text(size = text_size_legend, color = text_color),
      legend.text = element_text(size = text_size_legend, color = text_color),
      panel.grid = element_blank(),
      axis.title = element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "vertical"
    )
  
  return(p)
}
draw_custom_heatmap_ggplot_v3 <- function(data_list,
                                          data_type,
                                          p_thresholds = c(0.05, 0.01),
                                          color_rects = c("green", "red"),
                                          overlap_threshold = 3,
                                          apply_filling = TRUE,
                                          color_filling = "green",
                                          size_filling = 1,
                                          alpha_filling = 1,
                                          col_only_n_genes = FALSE,
                                          col_mapping_vector = NULL,
                                          row_mapping_vector = NULL,
                                          gene_list_sizes = TRUE,
                                          col_significant = FALSE,
                                          color_scale_range = NULL,
                                          palette = c("white", "#f1bcbb", "#edacab", "#e68a89"),
                                          text_color = "black",
                                          axis_text_angle = 45,
                                          axis_text_hjust = 0,
                                          text_size_tile = 3,
                                          text_size_axis = 10,
                                          text_size_legend = 10,
                                          title = NULL,
                                          title_size = 14) {
  
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  
  data <- data_list[[data_type]]$list
  cols_to_filter <- data_list[[data_type]]$cols
  rows_to_filter <- data_list[[data_type]]$rows
  
  if (col_significant) {
    data <- lapply(data, function(mat) mat[rows_to_filter, cols_to_filter])
  }
  
  chi2 <- data$chi2_value_matrix
  overlap <- data$number_overlap_matrix
  pval <- data$p_value_matrix
  
  if (gene_list_sizes) {
    if (is.null(col_mapping_vector)) {
      col_mapping_vector <- setNames(colnames(chi2), colnames(chi2))
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
      row_mapping_vector <- setNames(rownames(chi2), rownames(chi2))
    }
    row_mapping_vector <- setNames(
      paste0(row_mapping_vector, " (", data_list$gene_list_sizes[names(row_mapping_vector)], ")"),
      names(row_mapping_vector)
    )
  }
  
  df <- expand.grid(row = rownames(chi2), col = colnames(chi2), stringsAsFactors = FALSE) %>%
    mutate(
      chi2 = as.vector(chi2),
      overlap = as.vector(overlap),
      pval = as.vector(pval),
      row_label = row_mapping_vector[row],
      col_label = col_mapping_vector[col],
      label = sprintf("%.1f (%d)", chi2, overlap),
      log_chi2 = log2(chi2 + 1)
    )
  
  for (i in seq_along(p_thresholds)) {
    df[[paste0("sig_", i)]] <- df$pval < p_thresholds[i] & df$overlap >= overlap_threshold
  }
  
  if (is.null(color_scale_range)) {
    color_scale_range <- range(df$log_chi2, na.rm = TRUE)
  }
  
  p <- ggplot(df, aes(x = col_label, y = row_label)) +
    geom_tile(aes(fill = log_chi2), color = "white") +
    geom_text(aes(label = label), size = text_size_tile, color = text_color)
  
  for (i in seq_along(p_thresholds)) {
    sig_df <- df %>% filter(!!sym(paste0("sig_", i)))
    if (nrow(sig_df) > 0) {
      p <- p + geom_tile(data = sig_df, color = color_rects[i], fill = NA, size = 1.2)
    }
  }
  
  if (apply_filling) {
    p <- p + geom_point(data = df %>% filter(sig_1), shape = 16, size = size_filling,
                        color = color_filling, alpha = alpha_filling)
  }
  
  p <- p + geom_rect(data = NULL,
                     aes(xmin = 0.5, xmax = length(unique(df$col_label)) + 0.5,
                         ymin = 0.5, ymax = length(unique(df$row_label)) + 0.5),
                     inherit.aes = FALSE,
                     color = "black", fill = NA, size = 1)
  
  p <- p +
    scale_fill_gradientn(colors = palette,
                         limits = color_scale_range,
                         name = "log2(chi2 + 1)") +
    scale_x_discrete(position = "top") +
    scale_y_discrete(position = "right") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = axis_text_angle, hjust = axis_text_hjust,
                                 color = text_color, size = text_size_axis),
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
      theme(plot.title = element_text(hjust = 0, size = title_size, color = text_color))
  }
  
  return(p)
}
draw_custom_heatmap_ggplot_v3 <- function(data_list,
                                          data_type,
                                          p_thresholds = c(0.05, 0.01),
                                          color_rects = c("green", "red"),
                                          overlap_threshold = 3,
                                          apply_filling = TRUE,
                                          color_filling = "green",
                                          size_filling = 1,
                                          alpha_filling = 1,
                                          col_only_n_genes = FALSE,
                                          col_mapping_vector = NULL,
                                          row_mapping_vector = NULL,
                                          gene_list_sizes = TRUE,
                                          col_significant = FALSE,
                                          color_scale_range = NULL,
                                          palette = c("white", "#f1bcbb", "#edacab", "#e68a89"),
                                          text_color = "black",
                                          axis_text_angle = 45,
                                          axis_text_hjust = 0,
                                          text_size_tile = 3,
                                          text_size_axis = 10,
                                          text_size_legend = 10,
                                          title = NULL,
                                          title_size = 14) {
  
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  
  data <- data_list[[data_type]]$list
  cols_to_filter <- data_list[[data_type]]$cols
  rows_to_filter <- data_list[[data_type]]$rows
  
  if (col_significant) {
    data <- lapply(data, function(mat) mat[rows_to_filter, cols_to_filter])
  }
  
  chi2 <- data$chi2_value_matrix
  overlap <- data$number_overlap_matrix
  pval <- data$p_value_matrix
  
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
  
  df <- expand.grid(row = rownames(chi2), col = colnames(chi2), stringsAsFactors = FALSE) %>%
    mutate(
      chi2 = as.vector(chi2),
      overlap = as.vector(overlap),
      pval = as.vector(pval),
      row_label = row_mapping_vector[row],
      log_chi2 = log2(chi2 + 1),
      label = sprintf("%.1f (%d)", chi2, overlap)
    )
  
  sorted_col_order <- sort(cols_to_filter)
  df$col_factor <- factor(df$col, levels = sorted_col_order)
  
  for (i in seq_along(p_thresholds)) {
    df[[paste0("sig_", i)]] <- df$pval < p_thresholds[i] & df$overlap >= overlap_threshold
  }
  
  if (is.null(color_scale_range)) {
    color_scale_range <- range(df$log_chi2, na.rm = TRUE)
  }
  
  p <- ggplot(df, aes(x = col_factor, y = row_label)) +
    geom_tile(aes(fill = log_chi2), color = "white") +
    geom_text(aes(label = label), size = text_size_tile, color = text_color)
  
  for (i in seq_along(p_thresholds)) {
    sig_df <- df %>% filter(!!sym(paste0("sig_", i)))
    if (nrow(sig_df) > 0) {
      p <- p + geom_tile(data = sig_df, color = color_rects[i], fill = NA, size = 1.2)
    }
  }
  
  if (apply_filling) {
    p <- p + geom_point(data = df %>% filter(sig_1), shape = 16,
                        size = size_filling, color = color_filling, alpha = alpha_filling)
  }
  
  p <- p + geom_rect(data = NULL,
                     aes(xmin = 0.5, xmax = length(unique(df$col_factor)) + 0.5,
                         ymin = 0.5, ymax = length(unique(df$row_label)) + 0.5),
                     inherit.aes = FALSE, color = "black", fill = NA, size = 1)
  
  p <- p +
    scale_fill_gradientn(colors = palette,
                         limits = color_scale_range,
                         name = "log2(chi2 + 1)") +
    scale_x_discrete(position = "top",
                     labels = col_mapping_vector[sorted_col_order]) +
    scale_y_discrete(position = "right") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = axis_text_angle, hjust = axis_text_hjust,
                                 color = text_color, size = text_size_axis),
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
      theme(plot.title = element_text(hjust = 0, size = title_size, color = text_color))
  }
  
  return(p)
}

draw_custom_heatmap_ggplot_v3(
  data_list = grSignatures_factors_chi2Preprocessing,
  data_type = "significant_uniq_data",
  p_thresholds = c(0.05, 0.01),
  color_rects = c("#4C8D05", "#66023C"),
  overlap_threshold = 3,
  apply_filling = FALSE,
  col_only_n_genes = TRUE,
  col_significant = TRUE,
  color_scale_range = c(0, 4),
  palette = c("white", "#f8dedd", "#f1bcbb", "#edacab", "#e68a89"),
  
  # Nowe parametry:
  text_color = "black",            # kolor wszystkich tekstów
  axis_text_angle = 45,            # kąt etykiet osi X
  axis_text_hjust = 0,             # justowanie etykiet osi X
  text_size_tile = 3,              # wielkość tekstu w środku płytek
  text_size_axis = 10,             # wielkość etykiet osi
  text_size_legend = 10            # wielkość tekstu w legendzie
)
