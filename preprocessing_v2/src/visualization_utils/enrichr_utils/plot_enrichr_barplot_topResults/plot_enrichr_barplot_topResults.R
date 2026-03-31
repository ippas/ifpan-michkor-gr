plot_enrichr_barplot_topResults <- function(
    data,
    x_axis = "P.value",
    y_axis = "Term",
    n_genes_col = "n_genes",
    top_n = 15,
    title = NULL,
    x_label = expression(-log[10](P.value)),
    y_label = NULL,
    fill_color = "grey70",
    bar_width = 0.7,
    label_suffix = " genes",
    label_hjust = -0.3,
    title_size = 16,
    axis_label_size = 14,
    axis_text_size = 12,
    label_size_px = 16  # 👈 user input in pixels
) {
  library(dplyr)
  library(ggplot2)
  
  # --- 1️⃣ Validate input columns
  required_cols <- c(x_axis, y_axis)
  missing_cols <- required_cols[!required_cols %in% names(data)]
  if (length(missing_cols) > 0)
    stop(paste("Missing columns:", paste(missing_cols, collapse = ", ")))
  
  has_n_genes <- n_genes_col %in% names(data)
  
  # --- 2️⃣ Convert px → pt (1 px = 0.75 pt)
  label_size_pt <- label_size_px / 2.845
  
  # --- 3️⃣ Transform P/FDR values to -log10
  df <- data %>%
    mutate(log_value = -log10(.data[[x_axis]])) %>%
    arrange(.data[[x_axis]]) %>%
    slice_head(n = top_n)
  
  # --- 4️⃣ Base plot
  p <- ggplot(df, aes(x = log_value, y = reorder(.data[[y_axis]], log_value))) +
    geom_col(fill = fill_color, color = "black", width = bar_width)
  
  # --- 5️⃣ Add gene count labels (if available)
  if (has_n_genes) {
    p <- p + geom_text(
      aes(label = paste0(.data[[n_genes_col]], label_suffix)),
      hjust = label_hjust,
      size = label_size_pt,
      color = "black"
    )
  }
  
  # --- 6️⃣ Final styling
  p +
    labs(
      x = x_label,
      y = y_label,
      title = title
    ) +
    theme_classic(base_size = 14) +
    coord_cartesian(clip = "off") +
    theme(
      plot.title = element_text(size = title_size, hjust = 0.5, face = "bold"),
      axis.title.x = element_text(size = axis_label_size, face = "bold"),
      axis.title.y = element_text(size = axis_label_size, face = "bold"),
      axis.text.y = element_text(size = axis_text_size, color = "black"),
      axis.text.x = element_text(size = axis_text_size, color = "black"),
      plot.margin = margin(5, 70, 5, 5)
    )
}
