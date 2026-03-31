draw_custom_heatmap_fisher_ggplot_v3 <- function(
    data_list,
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
    palette = c("#2171b5", "#bdd7e7", "white", "#fcae91", "#cb181d"),
    text_color = "black",
    axis_text_angle = 45,
    axis_text_hjust = 0,
    text_size_tile = 3,
    text_size_axis = 10,
    text_size_legend = 10,
    title = NULL,
    title_size = 14
) {
  data <- data_list[[data_type]]$list
  cols_to_filter <- data_list[[data_type]]$cols
  rows_to_filter <- data_list[[data_type]]$rows
  
  if (col_significant) {
    data <- lapply(data, function(mat) mat[rows_to_filter, cols_to_filter])
  }
  
  fisher <- data$fisher_value_matrix
  overlap <- data$number_overlap_matrix
  pval <- data$p_value_matrix
  
  if (gene_list_sizes) {
    if (is.null(col_mapping_vector)) {
      col_mapping_vector <- stats::setNames(cols_to_filter, cols_to_filter)
    }
    col_mapping_vector <- stats::setNames(
      if (col_only_n_genes) {
        as.character(data_list$gene_list_sizes[names(col_mapping_vector)])
      } else {
        paste0(col_mapping_vector, " (", data_list$gene_list_sizes[names(col_mapping_vector)], ")")
      },
      names(col_mapping_vector)
    )
    
    if (is.null(row_mapping_vector)) {
      row_mapping_vector <- stats::setNames(rows_to_filter, rows_to_filter)
    }
    row_mapping_vector <- stats::setNames(
      paste0(row_mapping_vector, " (", data_list$gene_list_sizes[names(row_mapping_vector)], ")"),
      names(row_mapping_vector)
    )
  }
  
  df <- tidyr::expand_grid(
    row = rownames(fisher),
    col = colnames(fisher)
  ) %>%
    dplyr::mutate(
      fisher_value = as.vector(fisher),
      overlap = as.vector(overlap),
      pval = as.vector(pval),
      row_label = row_mapping_vector[row],
      col_factor = factor(col, levels = sort(cols_to_filter)),
      fisher_value = dplyr::if_else(is.na(fisher_value) | fisher_value <= 0, NA_real_, fisher_value),
      log_fisher = log2(fisher_value),
      label = dplyr::if_else(!is.na(log_fisher),
                             sprintf("%.1f (%d)", log_fisher, overlap),
                             "")
    )
  
  for (i in seq_along(p_thresholds)) {
    df[[paste0("sig_", i)]] <- df$pval < p_thresholds[i] & df$overlap >= overlap_threshold
  }
  
  # Automatyczny symetryczny zakres kolorów
  if (is.null(color_scale_range)) {
    finite_vals <- df$log_fisher[is.finite(df$log_fisher)]
    max_abs <- max(abs(finite_vals), na.rm = TRUE)
    color_scale_range <- c(-max_abs, max_abs)
  }
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = col_factor, y = row_label)) +
    ggplot2::geom_tile(ggplot2::aes(fill = log_fisher), color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = label), size = text_size_tile, color = text_color)
  
  for (i in seq_along(p_thresholds)) {
    sig_df <- dplyr::filter(df, !!rlang::sym(paste0("sig_", i)))
    if (nrow(sig_df) > 0) {
      p <- p + ggplot2::geom_tile(data = sig_df, color = color_rects[i], fill = NA, size = 1.2)
    }
  }
  
  if (apply_filling) {
    p <- p + ggplot2::geom_point(
      data = dplyr::filter(df, sig_1),
      shape = 16,
      size = size_filling,
      color = color_filling,
      alpha = alpha_filling
    )
  }
  
  p <- p + ggplot2::geom_rect(
    data = NULL,
    ggplot2::aes(
      xmin = 0.5, xmax = length(unique(df$col_factor)) + 0.5,
      ymin = 0.5, ymax = length(unique(df$row_label)) + 0.5
    ),
    inherit.aes = FALSE,
    color = "black", fill = NA, size = 1
  )
  
  p <- p +
    ggplot2::scale_fill_gradientn(
      colors = palette,
      values = scales::rescale(c(-2, 0, 2)),
      limits = color_scale_range,
      name = "log2(odds ratio)"
    ) +
    ggplot2::scale_x_discrete(
      position = "top",
      labels = col_mapping_vector[levels(df$col_factor)]
    ) +
    ggplot2::scale_y_discrete(position = "right") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = axis_text_angle,
        hjust = axis_text_hjust,
        color = text_color,
        size = text_size_axis
      ),
      axis.text.y = ggplot2::element_text(color = text_color, size = text_size_axis),
      legend.title = ggplot2::element_text(size = text_size_legend, color = text_color),
      legend.text = ggplot2::element_text(size = text_size_legend, color = text_color),
      panel.grid = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "vertical"
    )
  
  if (!is.null(title)) {
    p <- p + ggplot2::ggtitle(title) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          hjust = 0,
          size = title_size,
          color = text_color
        )
      )
  }
  
  return(p)
}
