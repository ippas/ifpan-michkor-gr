plot_icd10_bar_facet <- function(
    df,
    value_col,
    group_col = "signatureType",
    facet_col = "signatureType",
    fill_colors = NULL,
    group_order = NULL,
    order_within_group = "systemic",
    x_range = NULL,
    y_order = NULL,
    zero_one_normalize = FALSE,
    alpha = 1,
    border_color = NA,
    border_size = 0,
    n_cols = 2,
    facet_scales = "fixed",         # "fixed", "free_x", "free_y", "free"
    facet_strip_position = "top",   # "top", "bottom", "left", "right"
    show_legend = FALSE,
    
    # --- NEW: custom x-axis ticks/labels ---
    x_breaks = NULL,                # numeric vector, e.g. c(0,1,2,3)
    x_labels = NULL                 # labels vector or function, e.g. c("0","1","2","3")
) {
  
  value_col <- rlang::ensym(value_col)
  group_col <- rlang::ensym(group_col)
  facet_col <- rlang::ensym(facet_col)
  
  value_name <- rlang::as_name(value_col)
  group_name <- rlang::as_name(group_col)
  
  # --- ustal poziomy grup (priorytet: group_order > names(fill_colors) > dane) ---
  if (!is.null(group_order)) {
    group_levels <- as.character(group_order)
  } else if (!is.null(fill_colors) && !is.null(names(fill_colors)) && all(names(fill_colors) != "")) {
    group_levels <- as.character(names(fill_colors))
  } else {
    group_levels <- df %>%
      dplyr::pull(!!group_col) %>%
      unique() %>%
      as.character()
  }
  
  # --- factor dla group_col wg group_levels ---
  df <- df %>%
    dplyr::mutate(
      !!group_col := factor(as.character(!!group_col), levels = group_levels)
    )
  
  # --- opcjonalna normalizacja 0–1 w ramach grup ---
  if (zero_one_normalize) {
    df <- df %>%
      dplyr::group_by(!!group_col) %>%
      dplyr::mutate(
        !!value_col := {
          v <- .data[[value_name]]
          vmin <- min(v, na.rm = TRUE)
          vmax <- max(v, na.rm = TRUE)
          if (isTRUE(all.equal(vmin, vmax))) 0 else (v - vmin) / (vmax - vmin)
        }
      ) %>%
      dplyr::ungroup()
  }
  
  # --- domyślny y_order: na bazie wskazanej grupy ---
  if (is.null(y_order)) {
    order_levels <- df %>%
      dplyr::filter(as.character(.data[[group_name]]) == order_within_group) %>%
      dplyr::arrange(.data[[value_name]]) %>%
      dplyr::pull(mapped_icd10_range) %>%
      as.character()
  } else {
    order_levels <- as.character(y_order)
  }
  
  df_plot <- df %>%
    dplyr::mutate(
      mapped_icd10_range = factor(as.character(mapped_icd10_range), levels = order_levels)
    )
  
  # --- fill_colors: dopasowanie + uzupełnianie braków + wymuszenie kolejności wg group_levels ---
  if (is.null(fill_colors)) {
    fill_colors <- stats::setNames(
      grDevices::hcl.colors(length(group_levels), "Dark 3"),
      group_levels
    )
  } else {
    # jeśli brak nazw lub puste -> przypnij do group_levels po kolei
    if (is.null(names(fill_colors)) || any(names(fill_colors) == "")) {
      fill_colors <- stats::setNames(as.character(fill_colors), group_levels[seq_along(fill_colors)])
    }
    
    # uzupełnij brakujące poziomy
    missing_levels <- setdiff(group_levels, names(fill_colors))
    if (length(missing_levels) > 0) {
      extra_cols <- grDevices::hcl.colors(length(missing_levels), "Dark 3")
      names(extra_cols) <- missing_levels
      fill_colors <- c(fill_colors, extra_cols)
    }
    
    # finalne uporządkowanie
    fill_colors <- fill_colors[group_levels]
  }
  
  # --- wykres ---
  p <- ggplot2::ggplot(
    df_plot,
    ggplot2::aes(x = !!value_col, y = mapped_icd10_range, fill = !!group_col)
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = 0.7),
      width = 0.6,
      alpha = alpha,
      color = border_color,
      size = border_size
    ) +
    ggplot2::scale_fill_manual(values = fill_colors, drop = FALSE) +
    ggplot2::labs(
      x = rlang::as_label(value_col),
      y = "ICD-10 category",
      fill = rlang::as_label(group_col)
    ) +
    ggplot2::theme_classic() +
    theme_black_text +
    ggplot2::facet_wrap(
      ggplot2::vars(!!facet_col),
      ncol = n_cols,
      scales = facet_scales,
      strip.position = facet_strip_position
    )
  
  # --- x scale (with optional custom breaks/labels) ---
  if (!is.null(x_range)) {
    p <- p + ggplot2::scale_x_continuous(
      position = "top",
      limits = x_range,
      breaks = x_breaks,
      labels = x_labels
    )
  } else {
    p <- p + ggplot2::scale_x_continuous(
      position = "top",
      breaks = x_breaks,
      labels = x_labels
    )
  }
  
  if (!isTRUE(show_legend)) {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  
  return(p)
}


plot_icd10_bar_facet(
  preprocessing_summaryTable_grSystemicTissues_icd10F_df,
  fill_colors = c(
    systemic = "#335C67",
    neural   = "#540B0E",
    blood    = "#9E2A2B",
    lung     = "#E09F3E"
  ),   # przygaszona czerwień),
  mean_signif_cs,
  border_color = "black", n_cols = 4,
  border_size = 1.2,
  x_breaks = c(0, 2.5, 5),                # numeric vector, e.g. c(0,1,2,3)
  x_labels = c(0, 2.5, 5),  
  # x_range = c(0,5),
  y_order = rev(c("F3x", "F2x", "F0x", "F1x", "F4x", "F8x", "F5x", "F6x", "F7x", "F9x"))
) -> p1


plot_icd10_bar_facet(
  preprocessing_summaryTable_grSystemicTissues_icd10F_df,
  fill_colors = c(
    systemic = "#C9CBA3",
    neural   = "#FFE1A8",
    blood    = "#E26D5C",
    lung     = "#723D46"
  ),   # przygaszona czerwień),
  mean_signif_cs,
  border_color = "black", n_cols = 4,
  border_size = 1.2,
  x_breaks = c(0, 2.5, 5),                # numeric vector, e.g. c(0,1,2,3)
  x_labels = c(0, 2.5, 5),  
  # x_range = c(0,5),
  y_order = rev(c("F3x", "F2x", "F0x", "F1x", "F4x", "F8x", "F5x", "F6x", "F7x", "F9x"))
) -> p2


plot_icd10_bar_facet(
  preprocessing_summaryTable_grSystemicTissues_icd10F_df,
  fill_colors = c(
    systemic = "#335c67",
    neural   = "#fff3b0",
    blood    = "#9e2a2b",
    lung     = "#540b0e"
  ),   # przygaszona czerwień),
  mean_signif_cs,
  border_color = "black", n_cols = 4,
  border_size = 1.2,
  x_breaks = c(0, 2.5, 5),                # numeric vector, e.g. c(0,1,2,3)
  x_labels = c(0, 2.5, 5),  
  # x_range = c(0,5),
  y_order = rev(c("F3x", "F2x", "F0x", "F1x", "F4x", "F8x", "F5x", "F6x", "F7x", "F9x"))
) -> p3

plot_icd10_bar_facet(
  preprocessing_summaryTable_grSystemicTissues_icd10F_df,
  fill_colors = c(
    systemic = "#355070",
    neural   = "#b56576",
    blood    = "#e56b6f",
    lung     = "#eaac8b"
  ),   # przygaszona czerwień),
  mean_signif_cs,
  border_color = "black", n_cols = 4,
  border_size = 1.2,
  x_breaks = c(0, 2.5, 5),                # numeric vector, e.g. c(0,1,2,3)
  x_labels = c(0, 2.5, 5),  
  # x_range = c(0,5),
  y_order = rev(c("F3x", "F2x", "F0x", "F1x", "F4x", "F8x", "F5x", "F6x", "F7x", "F9x"))
) -> p4


plot_icd10_bar_facet(
  preprocessing_summaryTable_grSystemicTissues_icd10F_df,
  fill_colors = c(
    systemic = "#548687",
    neural   = "#ffffc7",
    blood    = "#b0413e",
    lung     = "#fcaa67"
  ),   # przygaszona czerwień),
  mean_signif_cs,
  border_color = "black", n_cols = 4,
  border_size = 1.2,
  x_breaks = c(0, 2.5, 5),                # numeric vector, e.g. c(0,1,2,3)
  x_labels = c(0, 2.5, 5),  
  # x_range = c(0,5),
  y_order = rev(c("F3x", "F2x", "F0x", "F1x", "F4x", "F8x", "F5x", "F6x", "F7x", "F9x"))
) -> p5

plot_icd10_bar_facet(
  preprocessing_summaryTable_grSystemicTissues_icd10F_df,
  fill_colors = c(
    systemic = "#0d1b2a",
    neural   = "#415a77",
    blood    = "#778da9",
    lung     = "#e0e1dd"
  ),   # przygaszona czerwień),
  mean_signif_cs,
  border_color = "black", n_cols = 4,
  border_size = 1.2,
  x_breaks = c(0, 2.5, 5),                # numeric vector, e.g. c(0,1,2,3)
  x_labels = c(0, 2.5, 5),  
  # x_range = c(0,5),
  y_order = rev(c("F3x", "F2x", "F0x", "F1x", "F4x", "F8x", "F5x", "F6x", "F7x", "F9x"))
) -> p6

p1 + p2 + p3 + p4 + p5 + p6
