# ##############################################################################
# ---- uses data ----
# ##############################################################################
flat_allGrSignatures_17.10.2025


# ##############################################################################
# ---- uses functions ----
# ##############################################################################
plot_full_overlap_heatmap

perform_chi2_tests # tutaj dodać -> dodane
processing_overlap_results


tmp <- run_full_overlap_analysis(
  gene_lists = flat_allGrSignatures_17.10.2025,
  total_genes = hgnc_symbols_vector_v110,
  plot_title = "log2(OR) overlap heatmap of GR-dependent gene sets",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  text_size_axis = 14,
  palette = c("navy", "white", "firebrick3")
)


tmp$processed$original_data$list$log2_odds_ratio_matrix

tmp$processed$original_data$df %>% filter(fdr < 0.0001) %>% 
  filter(!is.na(odds_ratio))



tmp_raw <- run_full_overlap_analysis(
  gene_lists = flat_allGrSignatures_17.10.2025[c("NeuralCellsUp",
                                                 "BloodCellsUp",
                                                 "LungCellsUp",
                                                 "NeuralCellsDown",
                                                 "BloodCellsDown",
                                                 "LungCellsDown",
                                                 "global_GR_genes_globalUp5TissuesDerivedCells",
                                                 "global_GR_genes_globalDown5TissuesDerivedCells")],
  total_genes = hgnc_symbols_vector_v110,
  plot_title = "",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  text_size_axis = 14,
  show_dendrograms = FALSE,
  scale_color_limit = c(-5, 5),
  text_contrast_range = c(-3, 4.9),
  p_thresholds = c(0.01, 0.0001),
  color_rects = c("#97C426", "#2F4603"),
  triangle_mode = "upper",
  palette = c("#07243e", "white", "darkred")
)

tmp2$plot



tmp2 <- run_full_overlap_analysis(
  gene_lists = flat_allGrSignatures_17.10.2025[c("minusGlobalUp5TissuesDerivedCells_NeuralCellsUp",
                                                 "minusGlobalUp5TissuesDerivedCells_BloodCellsUp",
                                                 "minusGlobalUp5TissuesDerivedCells_LungCellsUp",
                                                 "minusGlobalDown5TissuesDerivedCells_NeuralCellsDown",
                                                 "minusGlobalDown5TissuesDerivedCells_BloodCellsDown",
                                                 "minusGlobalDown5TissuesDerivedCells_LungCellsDown",
                                                 "global_GR_genes_globalUp5TissuesDerivedCells",
                                                 "global_GR_genes_globalDown5TissuesDerivedCells")],
  total_genes = hgnc_symbols_vector_v110,
  plot_title = "",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  text_size_axis = 14,
  show_dendrograms = FALSE,
  scale_color_limit = c(-5, 5),
  text_contrast_range = c(-3, 4.9),
  p_thresholds = c(0.01, 0.0001),
  color_rects = c("#97C426", "#2F4603"),
  triangle_mode = "upper",
  palette = c("#07243e", "white", "darkred")
)

tmp2$plot

tmp2 <- run_full_overlap_analysis(
  gene_lists = flat_allGrSignatures_17.10.2025[c("minusGlobalUp6TissuesDerivedCells_NeuralCellsUp",
                                                 "minusGlobalUp6TissuesDerivedCells_BloodCellsUp",
                                                 "minusGlobalUp6TissuesDerivedCells_LungCellsUp",
                                                 "minusGlobalDown6TissuesDerivedCells_NeuralCellsDown",
                                                 "minusGlobalDown6TissuesDerivedCells_BloodCellsDown",
                                                 "minusGlobalDown6TissuesDerivedCells_LungCellsDown",
                                                 "global_GR_genes_globalUp6TissuesDerivedCells",
                                                 "global_GR_genes_globalDown6TissuesDerivedCells")],
  total_genes = hgnc_symbols_vector_v110,
  plot_title = "",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  text_size_axis = 14,
  show_dendrograms = FALSE,
  scale_color_limit = c(-5, 5),
  text_contrast_range = c(-3, 4.9),
  p_thresholds = c(0.01, 0.0001),
  color_rects = c("#97C426", "#2F4603"),
  triangle_mode = "upper",
  palette = c("#07243e", "white", "darkred")
)

tmp2$plot


tmp2 <- run_full_overlap_analysis(
  gene_lists = flat_allGrSignatures_17.10.2025[c("minusGlobalUp4TissuesDerivedCells_NeuralCellsUp",
                                                 "minusGlobalUp4TissuesDerivedCells_BloodCellsUp",
                                                 "minusGlobalUp4TissuesDerivedCells_LungCellsUp",
                                                 "minusGlobalDown4TissuesDerivedCells_NeuralCellsDown",
                                                 "minusGlobalDown4TissuesDerivedCells_BloodCellsDown",
                                                 "minusGlobalDown4TissuesDerivedCells_LungCellsDown",
                                                 "global_GR_genes_globalUp4TissuesDerivedCells",
                                                 "global_GR_genes_globalDown4TissuesDerivedCells")],
  total_genes = hgnc_symbols_vector_v110,
  plot_title = "",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  text_size_axis = 14,
  show_dendrograms = FALSE,
  scale_color_limit = c(-5, 5),
  text_contrast_range = c(-3, 4.9),
  p_thresholds = c(0.01, 0.0001),
  color_rects = c("#97C426", "#2F4603"),
  triangle_mode = "upper",
  palette = c("#07243e", "white", "darkred")
)

tmp2$plot


heatmap_overlap_log2CHI2_ggplot(
  data_list = tmp2$processed,
  data_type = "original_data",
  # text_contrast_range = c(1, 3),
  # p_thresholds = c(0.05),
  # color_rects = "green",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  text_size_axis = 14,
  show_dendrograms = FALSE,
  # scale_color_limit = c(-5, 5),
  # text_contrast_range = c(0, 4.9),
  p_thresholds = c(0.01, 0.0001),
  color_rects = c("#97C426", "#2F4603"),
  triangle_mode = "upper",
  palette = c("white", "darkred")
)


tmp2$processed$original_data$df %>% 
  filter(!is.na(odds_ratio)) %>% 
  select(-c(fdr_value, fdr)) %>% 
  mutate(pair = pmap_chr(list(Var1, Var2), ~paste(sort(c(...)), collapse = "_vs_"))) %>%
  distinct(pair, .keep_all = TRUE) %>%
  select(-pair) %>% .$odds_ratio %>% boxplot


tmp_raw$processed$original_data$df %>% 
  filter(!is.na(odds_ratio)) %>% 
  select(-c(fdr_value, fdr)) %>% 
  mutate(pair = pmap_chr(list(Var1, Var2), ~paste(sort(c(...)), collapse = "_vs_"))) %>%
  distinct(pair, .keep_all = TRUE) %>%
  select(-pair) %>% .$odds_ratio %>% boxplot



bind_rows(
  tmp2$processed$original_data$df %>%
    filter(!is.na(odds_ratio)) %>%
    select(-c(fdr_value, fdr)) %>%
    mutate(
      pair = pmap_chr(list(Var1, Var2), ~paste(sort(c(...)), collapse = "_vs_")),
      source = "tmp2"
    ) %>%
    distinct(pair, .keep_all = TRUE),
  
  tmp_raw$processed$original_data$df %>%
    filter(!is.na(odds_ratio)) %>%
    select(-c(fdr_value, fdr)) %>%
    mutate(
      pair = pmap_chr(list(Var1, Var2), ~paste(sort(c(...)), collapse = "_vs_")),
      source = "tmp_raw"
    ) %>%
    distinct(pair, .keep_all = TRUE)
) %>%
  ggplot(aes(x = source, y = odds_ratio, fill = source)) +
  geom_boxplot(outlier.shape = 16, outlier.alpha = 0.4) +
  theme_minimal(base_size = 14) +
  labs(
    x = NULL,
    y = "Odds ratio",
    title = "Porównanie rozkładów wartości odds_ratio"
  ) +
  scale_fill_manual(values = c("tmp2" = "#4575b4", "tmp_raw" = "#d73027"))




df_combined <- bind_rows(
  tmp2$processed$original_data$df %>%
    filter(!is.na(odds_ratio)) %>%
    select(-c(fdr_value, fdr)) %>%
    mutate(
      pair = pmap_chr(list(Var1, Var2), ~paste(sort(c(...)), collapse = "_vs_")),
      source = "Minus global signatures (4 tissues)"
    ) %>%
    distinct(pair, .keep_all = TRUE),
  
  tmp_raw$processed$original_data$df %>%
    filter(!is.na(odds_ratio)) %>%
    select(-c(fdr_value, fdr)) %>%
    mutate(
      pair = pmap_chr(list(Var1, Var2), ~paste(sort(c(...)), collapse = "_vs_")),
      source = "Raw signatures"
    ) %>%
    distinct(pair, .keep_all = TRUE)
)

ggplot(df_combined, aes(x = source, y = odds_ratio)) +
  geom_boxplot(fill = "gray80", color = "black", outlier.shape = 16, outlier.alpha = 0.5) +
  geom_signif(
    comparisons = list(c("Raw signatures", "Minus global signatures (4 tissues)")),
    test = "wilcox.test",
    map_signif_level = TRUE,
    textsize = 4
  ) +
  theme_minimal(base_size = 14) +
  labs(
    x = NULL,
    y = "Odds ratio",
    title = "Comparison of odds ratio distributions"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(size = 12),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
  )
