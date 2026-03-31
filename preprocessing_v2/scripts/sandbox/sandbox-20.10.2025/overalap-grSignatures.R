# ============================================================
# 📘 Full overlap analyses for GR-dependent gene signatures
# ============================================================

library(dplyr)
library(purrr)
library(ggplot2)
library(ggsignif)

# --- 1️⃣ Analiza bazowa (raw signatures) ---
tmp_raw <- run_full_overlap_analysis(
  gene_lists = flat_allGrSignatures_17.10.2025[c(
    "NeuralCellsUp", "BloodCellsUp", "LungCellsUp",
    "NeuralCellsDown", "BloodCellsDown", "LungCellsDown",
    "global_GR_genes_globalUp5TissuesDerivedCells",
    "global_GR_genes_globalDown5TissuesDerivedCells"
  )],
  total_genes = hgnc_symbols_vector_v110,
  
  # 🔹 tytuły wykresów (osobne dla OR i chi2)
  plot_title_or = "",
  plot_title_chi2 = "",
  
  # 🔹 klastrowanie i wygląd osi
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  text_size_axis = 14,
  show_dendrograms = FALSE,
  
  # 🔹 ustawienia dla log2(OR)
  palette_or = c("#07243e", "white", "darkred"),
  color_scale_range_or = c(-5, 5),
  text_contrast_range_or = c(-3, 4.9),
  
  # 🔹 ustawienia dla log2(χ² + 1)
  palette_chi2 = c("white", "darkred"),
  color_scale_range_chi2 = c(0, 12),
  text_contrast_range_chi2 = c(0, 4),
  p_thresholds_chi2 = c(0.01, 0.0001),
  color_rects_chi2 = c("#97C426", "#2F4603"),
  
  # 🔹 inne
  triangle_mode = "upper"
)

tmp_raw$plot_chi2

# --- 2️⃣ Minus global signatures (4 tissues) ---
tmp_cluster <- run_full_overlap_analysis(
  gene_lists = flat_allGrSignatures_17.10.2025[c(
    "minusClustersKPO_NeuralCellsUp",
    "minusClustersKPO_BloodCellsUp",
    "minusClustersKPO_LungCellsUp",
    "minusClusterD_NeuralCellsDown",
    "minusClusterD_BloodCellsDown",
    "minusClusterD_LungCellsDown",
    "global_GR_genes_globalUp5TissuesDerivedCells",
    "global_GR_genes_globalDown5TissuesDerivedCells"
  )],
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

# --- 2️⃣ Minus global signatures (4 tissues) ---
tmp_4 <- run_full_overlap_analysis(
  gene_lists = flat_allGrSignatures_17.10.2025[c(
    "minusGlobalUp4TissuesDerivedCells_NeuralCellsUp",
    "minusGlobalUp4TissuesDerivedCells_BloodCellsUp",
    "minusGlobalUp4TissuesDerivedCells_LungCellsUp",
    "minusGlobalDown4TissuesDerivedCells_NeuralCellsDown",
    "minusGlobalDown4TissuesDerivedCells_BloodCellsDown",
    "minusGlobalDown4TissuesDerivedCells_LungCellsDown",
    "global_GR_genes_globalUp4TissuesDerivedCells",
    "global_GR_genes_globalDown4TissuesDerivedCells"
  )],
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

# --- 3️⃣ Minus global signatures (5 tissues) ---
tmp_5 <- run_full_overlap_analysis(
  gene_lists = flat_allGrSignatures_17.10.2025[c(
    "minusGlobalUp5TissuesDerivedCells_NeuralCellsUp",
    "minusGlobalUp5TissuesDerivedCells_BloodCellsUp",
    "minusGlobalUp5TissuesDerivedCells_LungCellsUp",
    "minusGlobalDown5TissuesDerivedCells_NeuralCellsDown",
    "minusGlobalDown5TissuesDerivedCells_BloodCellsDown",
    "minusGlobalDown5TissuesDerivedCells_LungCellsDown",
    "global_GR_genes_globalUp5TissuesDerivedCells",
    "global_GR_genes_globalDown5TissuesDerivedCells"
  )],
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

# --- 4️⃣ Minus global signatures (6 tissues) ---
tmp_6 <- run_full_overlap_analysis(
  gene_lists = flat_allGrSignatures_17.10.2025[c(
    "minusGlobalUp6TissuesDerivedCells_NeuralCellsUp",
    "minusGlobalUp6TissuesDerivedCells_BloodCellsUp",
    "minusGlobalUp6TissuesDerivedCells_LungCellsUp",
    "minusGlobalDown6TissuesDerivedCells_NeuralCellsDown",
    "minusGlobalDown6TissuesDerivedCells_BloodCellsDown",
    "minusGlobalDown6TissuesDerivedCells_LungCellsDown",
    "global_GR_genes_globalUp6TissuesDerivedCells",
    "global_GR_genes_globalDown6TissuesDerivedCells"
  )],
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



tmp_6$plot
# ============================================================
# 📊 Połączenie wyników w jedną ramkę
# ============================================================

extract_odds_ratios <- function(obj, label) {
  obj$processed$original_data$df %>%
    filter(!is.na(odds_ratio)) %>%
    select(-c(fdr_value, fdr)) %>%
    mutate(
      pair = pmap_chr(list(Var1, Var2), ~paste(sort(c(...)), collapse = "_vs_")),
      source = label
    ) %>%
    distinct(pair, .keep_all = TRUE)
}

df_combined <- bind_rows(
  extract_odds_ratios(tmp_raw, "Raw signatures"),
  extract_odds_ratios(tmp_4, "Minus global signatures (4 tissues)"),
  extract_odds_ratios(tmp_5, "Minus global signatures (5 tissues)"),
  extract_odds_ratios(tmp_6, "Minus global signatures (6 tissues)"),
  extract_odds_ratios(tmp_cluster, "Minus clusters KPOD, BMC")
)

# ============================================================
# 📈 Boxplot z testami statystycznymi (Wilcoxon)
# ============================================================

ggplot(df_combined, aes(x = source, y = odds_ratio)) +
  geom_boxplot(fill = "gray85", color = "black", outlier.shape = 16, outlier.alpha = 0.5) +
  geom_signif(
    comparisons = list(
      c("Raw signatures", "Minus global signatures (4 tissues)"),
      c("Raw signatures", "Minus global signatures (5 tissues)"),
      c("Raw signatures", "Minus global signatures (6 tissues)"),
      c("Raw signatures", "Minus clusters KPOD, BMC")
    ),
    test = "wilcox.test",
    map_signif_level = TRUE,
    step_increase = 0.12,
    textsize = 4
  ) +
  theme_minimal(base_size = 14) +
  labs(
    x = NULL,
    y = "Odds ratio",
    title = "Comparison of odds ratio distributions across signature types"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(size = 12, angle = 30, hjust = 1),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
  )
