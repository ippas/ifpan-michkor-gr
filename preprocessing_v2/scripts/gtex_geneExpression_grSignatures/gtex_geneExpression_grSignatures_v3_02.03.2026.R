
hgnc_symbols_vector_v110_gtex %>% 
  select(ontologyId, geneSymbol, tissue, median_max)

hgnc_symbols_vector_v110_gtex %>%
  select(ontologyId, geneSymbol, tissue, median_max) %>%
  ggplot(aes(x = tissue, y = log2(median_max))) +
  geom_boxplot() +
  theme_classic() +
  labs(
    x = "Tissue",
    y = "Median max expression",
    title = "Distribution of median_max across tissues"
  )


df_plot <- hgnc_symbols_vector_v110_gtex %>%
  select(ontologyId, geneSymbol, tissue, median_max) %>%
  filter(!is.na(tissue)) %>%
  mutate(
    tissue = factor(tissue),
    tissue_id = as.numeric(tissue),
    stripe = tissue_id %% 2   # 0 / 1 naprzemiennie
  )


ggplot(df_plot, aes(x = tissue, y = log2(median_max))) +
  
  # boxplot bez outlierów
  geom_boxplot(outlier.shape = NA, width = 0.6, color = "black") +
  
  # jitter z naprzemiennym kolorem
  # geom_jitter(
  #   aes(color = factor(stripe)),
  #   width = 0.2,
  #   alpha = 0.1,
  #   size = 0.6
  # ) +
  
  scale_color_manual(
    values = c("0" = "#1f78b4", "1" = "#e31a1c"),
    guide = "none"
  ) +
  
  theme_classic() +
  labs(
    x = "Tissue",
    y = "log2(median max expression)"
  )


# ##############################################################################
plot_gene_3metrics <- function(gene_symbol,
                               tissues = NULL,
                               log2 = TRUE,
                               decreasing = TRUE,
                               data = hgnc_symbols_vector_v110_gtex) {
  
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  
  # --- filtrowanie ---
  df <- data %>%
    select(tissue, geneSymbol, median_max) %>%
    filter(!is.na(tissue))
  
  if (!is.null(tissues)) {
    df <- df %>% filter(tissue %in% tissues)
  }
  
  # --- RAW (log2) ---
  df_gene <- df %>%
    filter(geneSymbol == gene_symbol) %>%
    group_by(tissue) %>%
    summarise(raw_log2 = mean(median_max, na.rm = TRUE), .groups = "drop")
  
  if (nrow(df_gene) == 0) stop("Gene not found.")
  
  # --- tło ---
  bg <- df %>%
    group_by(tissue) %>%
    summarise(bg_median_log2 = median(median_max, na.rm = TRUE),
              .groups = "drop")
  
  # --- percentyl ---
  pct <- df %>%
    group_by(tissue) %>%
    mutate(percentile = percent_rank(median_max)) %>%
    ungroup() %>%
    filter(geneSymbol == gene_symbol) %>%
    group_by(tissue) %>%
    summarise(percentile = mean(percentile, na.rm = TRUE),
              .groups = "drop")
  
  # --- łączenie ---
  df_all <- df_gene %>%
    left_join(bg, by = "tissue") %>%
    left_join(pct, by = "tissue") %>%
    mutate(
      norm_vs_median = if (log2) {
        raw_log2 - bg_median_log2
      } else {
        2^(raw_log2 - bg_median_log2)
      }
    )
  
  # --- sortowanie ---
  if (decreasing) {
    df_all <- df_all %>% arrange(desc(norm_vs_median))
  } else {
    df_all <- df_all %>% arrange(norm_vs_median)
  }
  
  df_all <- df_all %>%
    mutate(tissue = factor(tissue, levels = tissue))
  
  # --- long ---
  df_long <- df_all %>%
    select(tissue, norm_vs_median, percentile, raw_log2) %>%
    pivot_longer(-tissue, names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(
      metric,
      levels = c("norm_vs_median", "percentile", "raw_log2"),
      labels = c(
        if (log2) "log2FC vs tissue median"
        else "Fold-change vs tissue median",
        "Percentile (0–1)",
        "Raw (log2)"
      )
    ))
  
  ggplot(df_long, aes(x = tissue, y = value)) +
    geom_col(width = 0.8) +
    coord_flip() +
    facet_wrap(~ metric, scales = "free_x", ncol = 3) +
    theme_classic() +
    labs(
      title = paste0(gene_symbol, " across tissues"),
      x = "Tissue",
      y = NULL
    )
}
plot_gene_3metrics("IL6", log2 = T, decreasing = F)


# ##############################################################################
plot_geneset_gtex_box <- function(
    genes,
    tissues = NULL,
    metrics = c("norm_vs_median", "percentile", "raw_log2"),
    log2 = TRUE,
    sort_by = "norm_vs_median",
    decreasing = TRUE,
    keep_original_order = FALSE,
    tissue_order = NULL,
    ncol = 3,
    coord_flip = TRUE,
    x_angle = 45,
    add.dots = FALSE,
    dot.alpha = 0.6,
    dot.size = 1.5,
    dot.width = 0.2,
    expression_limit = NULL,
    data = hgnc_symbols_vector_v110_gtex
) {
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  if (length(genes) == 0) stop("genes vector is empty.")
  
  df <- data %>%
    select(tissue, geneSymbol, median_max) %>%
    filter(!is.na(tissue))
  
  if (!is.null(tissues)) {
    df <- df %>% filter(tissue %in% tissues)
  }
  
  bg <- df %>%
    group_by(tissue) %>%
    summarise(bg_median_log2 = median(median_max, na.rm = TRUE),
              .groups = "drop")
  
  df_pct <- df %>%
    group_by(tissue) %>%
    mutate(percentile = percent_rank(median_max)) %>%
    ungroup()
  
  df_set <- df_pct %>%
    filter(geneSymbol %in% genes)
  
  if (nrow(df_set) == 0) stop("None of the provided genes were found.")
  
  df_set2 <- df_set %>%
    left_join(bg, by = "tissue") %>%
    mutate(
      raw_log2 = median_max,
      norm_vs_median = if (log2) {
        raw_log2 - bg_median_log2
      } else {
        2^(raw_log2 - bg_median_log2)
      }
    )
  
  df_long <- df_set2 %>%
    select(tissue, geneSymbol, raw_log2, norm_vs_median, percentile) %>%
    pivot_longer(
      cols = c(raw_log2, norm_vs_median, percentile),
      names_to = "metric",
      values_to = "value"
    ) %>%
    filter(metric %in% metrics)
  
  # --- ordering ---
  if (!is.null(tissue_order)) {
    df_long$tissue <- factor(df_long$tissue, levels = tissue_order)
  } else if (!keep_original_order) {
    ord <- df_long %>%
      filter(metric == sort_by) %>%
      group_by(tissue) %>%
      summarise(med = median(value, na.rm = TRUE), .groups = "drop") %>%
      arrange(if (decreasing) desc(med) else med) %>%
      pull(tissue)
    df_long$tissue <- factor(df_long$tissue, levels = ord)
  }
  
  p <- ggplot(df_long, aes(x = tissue, y = value)) +
    geom_boxplot(outlier.shape = NA, width = 0.7)
  
  if (add.dots) {
    p <- p +
      geom_jitter(width = dot.width,
                  alpha = dot.alpha,
                  size = dot.size)
  }
  
  if (length(metrics) > 1) {
    p <- p + facet_wrap(~ metric, scales = "free_y", ncol = ncol)
  }
  
  p <- p +
    theme_classic() +
    labs(
      title = paste0("Gene set across tissues (n=", 
                     length(unique(df_set2$geneSymbol)), ")"),
      x = "Tissue",
      y = NULL
    )
  
  if (!is.null(expression_limit)) {
    if (length(expression_limit) == 1) {
      p <- p + coord_cartesian(ylim = c(-expression_limit, expression_limit))
    }
    if (length(expression_limit) == 2) {
      p <- p + coord_cartesian(ylim = expression_limit)
    }
  }
  
  if (coord_flip) {
    p <- p + coord_flip()
  } else {
    p <- p +
      theme(axis.text.x = element_text(angle = x_angle, hjust = 1))
  }
  
  return(p)
}


plot_geneset_gtex_box(flat_allGrSignatures_31.10.2025$NeuralCellsUp, 
                      tissues = c("Lung", "Brain", "Whole_Blood"), log2 = T,
                      add.dots = T,
                      metrics = c("norm_vs_median", "raw_log2"),
                      # expression_limit = c(0, 300),
                      expression_limit = c(0, 300),
                      coord_flip = F)
