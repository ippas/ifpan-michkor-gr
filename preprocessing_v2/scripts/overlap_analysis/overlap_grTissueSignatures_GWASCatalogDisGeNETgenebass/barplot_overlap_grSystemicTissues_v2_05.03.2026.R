preprocessing_summaryTable_grSystemicTissues_icd10F_df %>% 
  group_by(mapped_icd10_range, tissue) %>% 
  nest %>%
  mutate(data_signif = map(data, ~ .x %>% filter(p_value < 0.05 & observed_overlap >= 3))) %>% 
  mutate(data_ns     = map(data, ~ .x %>% filter(p_value > 0.05 | observed_overlap < 3))) %>% 
  mutate(signif_cs   = map(data_signif, ~ .x$combine_score)) %>% 
  mutate(ns_cs       = map(data_ns, ~ .x$combine_score)) %>% 
  mutate(n_phenotypes          = map(data, ~ .x$phenotype %>% unique %>% length)) %>% unnest(n_phenotypes) %>% 
  mutate(n_signif_phenotypes   = map(data_signif, ~ .x$phenotype %>% unique %>% length)) %>% unnest(n_signif_phenotypes) %>% 
  mutate(n_signif_associations = map(data_signif, ~ .x %>% nrow)) %>% unnest(n_signif_associations) %>% 
  mutate(mean_signif_cs         = map(data_signif, ~ .x$combine_score %>% mean)) %>% 
  mutate(mean_all_cs            = map(data, ~ .x$combine_score %>% mean)) %>% 
  mutate(sum_signif_cs          = map(data_signif, ~ .x$combine_score %>% sum)) %>% 
  mutate(sum_ns_cs              = map(data_ns, ~ .x$combine_score %>% sum)) %>% 
  unnest(c(mean_signif_cs, mean_all_cs, sum_signif_cs, sum_ns_cs)) %>% 
  mutate(mean_signif_cs         = replace(mean_signif_cs, is.nan(mean_signif_cs), 0)) %>% 
  mutate(data_top3systemic      = map(data, ~ .x %>% filter(phenotype %in% top3grSystemic_phenotypes))) %>% 
  mutate(data_top3Tissues       = map(data, ~ .x %>% filter(phenotype %in% top3grTissues_phenotypes))) %>% 
  mutate(top3Tissues_sumSignif_cs = map(data_top3Tissues, ~.x %>% filter(p_value < 0.05 & observed_overlap >= 3) %>% .$combine_score %>% sum)) %>% 
  mutate(top3Tissues_sumNS_cs   = map(data_top3Tissues, ~ .x %>% filter(p_value > 0.05 | observed_overlap < 3) %>% .$combine_score %>% sum)) %>% 
  mutate(top3Systemic_sumSignif_cs    = map(data_top3systemic, ~.x %>% filter(p_value < 0.05 & observed_overlap >= 3) %>% .$combine_score %>% sum)) %>% 
  mutate(top3Systemic_sumNS_cs  = map(data_top3systemic, ~.x %>% filter(p_value > 0.05 | observed_overlap < 3) %>% .$combine_score %>% sum)) %>% 
  unnest(top3Systemic_sumSignif_cs, top3Systemic_sumNS_cs, top3Tissues_sumSignif_cs, top3Tissues_sumNS_cs) -> preprocessing_summaryTable_grSystemicTissues_icd10F_df

preprocessing_summaryTable_grSystemicTissues_icd10F_df %>% 
  rename(signatureType = "tissue") -> preprocessing_summaryTable_grSystemicTissues_icd10F_df

plot_icd10_bar <- function(
    df,
    value_col,
    group_col = "signatureType",
    fill_colors = NULL,
    group_order = NULL,
    order_within_group = "systemic",
    x_range = NULL,
    y_order = NULL,
    zero_one_normalize = FALSE,
    alpha = 1,
    border_color = NA,
    border_size = 0
) {
  
  value_col <- rlang::ensym(value_col)
  group_col <- rlang::ensym(group_col)
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
          if (isTRUE(all.equal(vmin, vmax))) {
            0
          } else {
            (v - vmin) / (vmax - vmin)
          }
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
    } else {
      fill_colors <- fill_colors
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
      position = ggplot2::position_dodge(width = 0.85),
      width = 0.8,
      alpha = alpha,
      color = border_color,
      linewidth = border_size
    ) +
    ggplot2::scale_fill_manual(values = fill_colors, drop = FALSE) +
    ggplot2::labs(
      x = rlang::as_label(value_col),
      y = "ICD-10 category",
      fill = rlang::as_label(group_col)
    ) +
    ggplot2::theme_classic() +
    theme_black_text
  
  if (!is.null(x_range)) {
    p <- p + ggplot2::scale_x_continuous(position = "top", limits = x_range)
  } else {
    p <- p + ggplot2::scale_x_continuous(position = "top")
  }
  
  return(p)
}


plot_icd10_bar(
  preprocessing_summaryTable_grSystemicTissues_icd10F_df %>% 
    mutate(mean_signif_cs_2 = ifelse(mean_signif_cs == 0, mean_signif_cs + 0.002, mean_signif_cs)),
  fill_colors = c(
    lung     = "#fcaa67",  
    blood    = "#b0413e",
    neural   = "#ffffc7",
    systemic = "#548687"),   # przygaszona czerwień),
  mean_signif_cs_2,
  border_color = NA,
  # x_range = c(0,5),
  y_order = rev(c("F3x", "F2x", "F0x", "F1x", "F4x", "F8x", "F5x", "F6x", "F7x", "F9x"))
) -> p1

plot_icd10_bar(
  preprocessing_summaryTable_grSystemicTissues_icd10F_df,
  fill_colors = c(
    lung     = "#E09F3E",  
    blood    = "#9E2A2B",
    neural   = "#540B0E",
    systemic = "#335C67"),   # przygaszona czerwień),
  mean_signif_cs,
  border_color = "black",
  # x_range = c(0,5),
  y_order = rev(c("F3x", "F2x", "F0x", "F1x", "F4x", "F8x", "F5x", "F6x", "F7x", "F9x"))
) -> p2



plot_icd10_bar(
  preprocessing_summaryTable_grSystemicTissues_icd10F_df,
  fill_colors = c(
    lung     = "#e0e1dd",  # stonowany niebieski
    blood    = "#778da9",
    neural   = "#415a77",
    systemic = "#0d1b2a"),   # przygaszona czerwień),
  mean_signif_cs,
  border_color = "black",
  # x_range = c(0,5),
  y_order = rev(c("F3x", "F2x", "F0x", "F1x", "F4x", "F8x", "F5x", "F6x", "F7x", "F9x"))
) -> p3

plot_icd10_bar(
  preprocessing_summaryTable_grSystemicTissues_icd10F_df,
  fill_colors = c(
    lung     = "#FFE1A8",  # stonowany niebieski
    blood    = "#E26D5C",
    neural   = "#723D46",
    systemic = "#C9CBA3"),   # przygaszona czerwień),
  mean_signif_cs,
  border_color = "black",
  # x_range = c(0,5),
  y_order = rev(c("F3x", "F2x", "F0x", "F1x", "F4x", "F8x", "F5x", "F6x", "F7x", "F9x"))
) -> p4



plot_icd10_bar(
  preprocessing_summaryTable_grSystemicTissues_icd10F_df,
  fill_colors = c(
    lung     = "#eaac8b",  # stonowany niebieski
    blood    = "#e56b6f",
    neural   = "#b56576",
    systemic = "#355070"),   # przygaszona czerwień),
  mean_signif_cs,
  border_color = "black",
  # x_range = c(0,5),
  y_order = rev(c("F3x", "F2x", "F0x", "F1x", "F4x", "F8x", "F5x", "F6x", "F7x", "F9x"))
) -> p5

plot_icd10_bar(
  preprocessing_summaryTable_grSystemicTissues_icd10F_df,
  fill_colors = c(
    lung     = "#d8572a",  # stonowany niebieski
    blood    = "#db7c26",
    neural   = "#f7b538",
    systemic = "#780116"),   # przygaszona czerwień),
  mean_signif_cs,
  border_color = "black",
  # x_range = c(0,5),
  y_order = rev(c("F3x", "F2x", "F0x", "F1x", "F4x", "F8x", "F5x", "F6x", "F7x", "F9x"))
) -> p6



p2 + p4  + p1 + p3 +  p5 + p6
 

svg("/home/mateusz/projects/ifpan-michkor-gr/results_v2/overlap/overlap_grTissues_GWASCatalogDisGeNETgenebass/figures/barplot_grTissuesSystemic_perTissue_meanCS_v2_04.03.2026.svg",
    width = 8,
    height = 8)

(p1 + p2 + p1 + p2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

dev.off()



