library(dplyr)
library(purrr)
library(ggplot2)
library(ggsignif)

# ============================================================
# 1️⃣ Raw GR-dependent signatures (baseline)
# ============================================================
overlap_raw_signatures <- run_full_overlap_analysis(
  gene_lists = flat_allGrSignatures_17.10.2025[c(
    "NeuralCellsUp", "BloodCellsUp", "LungCellsUp",
    "NeuralCellsDown", "BloodCellsDown", "LungCellsDown"
  )],
  total_genes = hgnc_symbols_vector_v110,
  
  plot_title_or = "       Raw GR-dependent signatures — log2(OR)",
  plot_title_chi2 = "       Raw GR-dependent signatures — log2(χ² + 1)",
  
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  text_size_axis = 14,
  show_dendrograms = FALSE,
  
  palette_or = c("#07243e", "white", "darkred"),
  color_scale_range_or = c(-5, 5),
  text_contrast_range_or = c(-3, 4.9),
  
  palette_chi2 = c("white", "darkred"),
  color_scale_range_chi2 = c(0, 8),
  text_contrast_range_chi2 = c(0, 6),
  p_thresholds_chi2 = c(0.01, 0.0001),
  color_rects_chi2 = c("#97C426", "#2F4603"),
  
  triangle_mode = "upper"
)

# ============================================================
# 2️⃣ Minus clusters (KPOD/BMC removed)
# ============================================================
# Mapper – usuwa prefiksy "minusClustersKPO_" i "minusClusterD_"
label_mapper_clusters <- c(
  "minusClustersKPO_NeuralCellsUp" = "NeuralCellsUp",
  "minusClustersKPO_BloodCellsUp"  = "BloodCellsUp",
  "minusClustersKPO_LungCellsUp"   = "LungCellsUp",
  "minusClusterD_NeuralCellsDown"  = "NeuralCellsDown",
  "minusClusterD_BloodCellsDown"   = "BloodCellsDown",
  "minusClusterD_LungCellsDown"    = "LungCellsDown"
)

overlap_minus_clusters <- run_full_overlap_analysis(
  gene_lists = flat_allGrSignatures_17.10.2025[c(
    "minusClustersKPO_NeuralCellsUp",
    "minusClustersKPO_BloodCellsUp",
    "minusClustersKPO_LungCellsUp",
    "minusClusterD_NeuralCellsDown",
    "minusClusterD_BloodCellsDown",
    "minusClusterD_LungCellsDown"
  )],
  total_genes = hgnc_symbols_vector_v110,
  
  plot_title_or = "Minus cluster-derived signatures — log2(OR)",
  plot_title_chi2 = "Minus cluster-derived signatures — log2(χ² + 1)",
  
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  text_size_axis = 14,
  show_dendrograms = FALSE,
  
  palette_or = c("#07243e", "white", "darkred"),
  color_scale_range_or = c(-5, 5),
  text_contrast_range_or = c(-3, 4.9),
  
  palette_chi2 = c("white", "darkred"),
  color_scale_range_chi2 = c(0, 8),
  text_contrast_range_chi2 = c(0, 6),
  p_thresholds_chi2 = c(0.01, 0.0001),
  color_rects_chi2 = c("#97C426", "#2F4603"),
  
  triangle_mode = "upper",
  
  # 👇 przekazanie mappera do obu heatmap
  row_labels_map = label_mapper_clusters,
  col_labels_map = label_mapper_clusters
)

# ============================================================
# 4⃣ Minus global (5 tissues)
# ============================================================

# Mapper — usuwa prefiksy "minusGlobalUp6TissuesDerivedCells_" i "minusGlobalDown6TissuesDerivedCells_"
label_mapper_5tissue <- c(
  "minusGlobalUp5TissuesDerivedCells_NeuralCellsUp"     = "NeuralCellsUp",
  "minusGlobalUp5TissuesDerivedCells_BloodCellsUp"      = "BloodCellsUp",
  "minusGlobalUp5TissuesDerivedCells_LungCellsUp"       = "LungCellsUp",
  "minusGlobalDown5TissuesDerivedCells_NeuralCellsDown" = "NeuralCellsDown",
  "minusGlobalDown5TissuesDerivedCells_BloodCellsDown"  = "BloodCellsDown",
  "minusGlobalDown5TissuesDerivedCells_LungCellsDown"   = "LungCellsDown"
)

# Analiza: minus global (6 tissues)
overlap_minus_global5 <- run_full_overlap_analysis(
  gene_lists = flat_allGrSignatures_17.10.2025[c(
    "minusGlobalUp5TissuesDerivedCells_NeuralCellsUp",
    "minusGlobalUp5TissuesDerivedCells_BloodCellsUp",
    "minusGlobalUp5TissuesDerivedCells_LungCellsUp",
    "minusGlobalDown5TissuesDerivedCells_NeuralCellsDown",
    "minusGlobalDown5TissuesDerivedCells_BloodCellsDown",
    "minusGlobalDown5TissuesDerivedCells_LungCellsDown"
  )],
  total_genes = hgnc_symbols_vector_v110,
  
  plot_title_or  = "Minus global (5-tissue criterion) — log2(OR)",
  plot_title_chi2 = "Minus global (5-tissue criterion) — log2(χ² + 1)",
  
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  text_size_axis = 14,
  show_dendrograms = FALSE,
  
  palette_or = c("#07243e", "white", "darkred"),
  color_scale_range_or = c(-5, 5),
  text_contrast_range_or = c(-3, 4.9),
  
  palette_chi2 = c("white", "darkred"),
  color_scale_range_chi2 = c(0, 8),
  text_contrast_range_chi2 = c(0, 5.8),
  p_thresholds_chi2 = c(0.01, 0.0001),
  color_rects_chi2 = c("#97C426", "#2F4603"),
  
  triangle_mode = "upper",
  
  # 👇 przekazanie mappera do obu heatmap
  row_labels_map = label_mapper_5tissue,
  col_labels_map = label_mapper_5tissue
)




# ============================================================
# 4️⃣ Minus global (6 tissues)
# ============================================================

# Mapper — usuwa prefiksy "minusGlobalUp6TissuesDerivedCells_" i "minusGlobalDown6TissuesDerivedCells_"
label_mapper_6tissue <- c(
  "minusGlobalUp6TissuesDerivedCells_NeuralCellsUp"     = "NeuralCellsUp",
  "minusGlobalUp6TissuesDerivedCells_BloodCellsUp"      = "BloodCellsUp",
  "minusGlobalUp6TissuesDerivedCells_LungCellsUp"       = "LungCellsUp",
  "minusGlobalDown6TissuesDerivedCells_NeuralCellsDown" = "NeuralCellsDown",
  "minusGlobalDown6TissuesDerivedCells_BloodCellsDown"  = "BloodCellsDown",
  "minusGlobalDown6TissuesDerivedCells_LungCellsDown"   = "LungCellsDown"
)

# Analiza: minus global (6 tissues)
overlap_minus_global6 <- run_full_overlap_analysis(
  gene_lists = flat_allGrSignatures_17.10.2025[c(
    "minusGlobalUp6TissuesDerivedCells_NeuralCellsUp",
    "minusGlobalUp6TissuesDerivedCells_BloodCellsUp",
    "minusGlobalUp6TissuesDerivedCells_LungCellsUp",
    "minusGlobalDown6TissuesDerivedCells_NeuralCellsDown",
    "minusGlobalDown6TissuesDerivedCells_BloodCellsDown",
    "minusGlobalDown6TissuesDerivedCells_LungCellsDown"
  )],
  total_genes = hgnc_symbols_vector_v110,
  
  plot_title_or  = "       Minus global (6-tissue criterion) — log2(OR)",
  plot_title_chi2 = "       Minus global (6-tissue criterion) — log2(χ² + 1)",
  
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  text_size_axis = 14,
  show_dendrograms = FALSE,
  
  palette_or = c("#07243e", "white", "darkred"),
  color_scale_range_or = c(-5, 5),
  text_contrast_range_or = c(-3, 4.9),
  
  palette_chi2 = c("white", "darkred"),
  color_scale_range_chi2 = c(0, 8),
  text_contrast_range_chi2 = c(0, 5.8),
  p_thresholds_chi2 = c(0.01, 0.0001),
  color_rects_chi2 = c("#97C426", "#2F4603"),
  
  triangle_mode = "upper",
  
  # 👇 przekazanie mappera do obu heatmap
  row_labels_map = label_mapper_6tissue,
  col_labels_map = label_mapper_6tissue
)


overlap_minus_global6$processed$original_data$df



# ##############################################################################
# ---- chi2 heatmap ----
# ##############################################################################
combined_plot <- wrap_plots(
  overlap_raw_signatures$plot_chi2,
  overlap_minus_clusters$plot_chi2,
  overlap_minus_global6$plot_chi2,
  overlap_minus_global5$plot_chi2,
  ncol = 2,
  guides = "collect"
) & theme(legend.position = "bottom")

# ============================================================
# 🟩🟥 3️⃣ Utworzenie niestandardowej legendy (puste prostokąty)
# ============================================================

custom_legend <- create_customRect_patch_legend(
  colors = c("#97C426", "#2F4603"),
  labels = c("p < 0.01", "p < 0.0001"),
  spacing = 2,
  box_linewidth = 4
)

# ============================================================
# 🧩 4️⃣ Połączenie wszystkiego (wykresy + legenda)
# ============================================================

final_plot <- wrap_plots(
  combined_plot,
  custom_legend,
  ncol = 1,
  heights = c(10, 1)
)

# ============================================================
# 💾 5️⃣ Podgląd lub zapis
# ============================================================

final_plot

dev.off()

# ##############################################################################
# ---- OR heatmap ----
# ##############################################################################
combined_plot <- wrap_plots(
  overlap_raw_signatures$plot_or,
  overlap_minus_clusters$plot_or,
  overlap_minus_global6$plot_or,
  overlap_minus_global5$plot_or,
  ncol = 2,
  guides = "collect"
) & theme(legend.position = "bottom")

# ============================================================
# 🟩🟥 3️⃣ Utworzenie niestandardowej legendy (puste prostokąty)
# ============================================================

custom_legend <- create_customRect_patch_legend(
  colors = c("#97C426", "#2F4603"),
  labels = c("p < 0.01", "p < 0.0001"),
  spacing = 2,
  box_linewidth = 4
)

# ============================================================
# 🧩 4️⃣ Połączenie wszystkiego (wykresy + legenda)
# ============================================================

final_plot <- wrap_plots(
  combined_plot,
  custom_legend,
  ncol = 1,
  heights = c(10, 1)
)

# ============================================================
# 💾 5️⃣ Podgląd lub zapis
# ============================================================

final_plot



# ##############################################################################
# ---- save to files ----
# ##############################################################################
# ============================================================
# 📁 Ścieżka zapisu wyników
# ============================================================
output_dir <- "/home/mateusz/projects/ifpan-michkor-gr/results_v2/positive_validation_grSignatures_21.10.2025/overlap_chi2"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ============================================================
# 💾 Zapis pliku SVG dla χ² heatmap
# ============================================================
combined_plot <- wrap_plots(
  overlap_raw_signatures$plot_chi2,
  overlap_minus_clusters$plot_chi2,
  overlap_minus_global6$plot_chi2,
  overlap_minus_global5$plot_chi2,
  ncol = 2,
  guides = "collect"
) & theme(legend.position = "bottom")

# ============================================================
# 🟩🟥 3️⃣ Utworzenie niestandardowej legendy (puste prostokąty)
# ============================================================

custom_legend <- create_customRect_patch_legend(
  colors = c("#97C426", "#2F4603"),
  labels = c("p < 0.01", "p < 0.0001"),
  spacing = 2,
  box_linewidth = 4
)

# ============================================================
# 🧩 4️⃣ Połączenie wszystkiego (wykresy + legenda)
# ============================================================

final_plot <- wrap_plots(
  combined_plot,
  custom_legend,
  ncol = 1,
  heights = c(10, 1)
)

# ============================================================
# 💾 5️⃣ Podgląd lub zapis
# ============================================================

final_plot
chi2_svg_path <- file.path(output_dir, "combined_overlap_chi2_finalPlot_22.10.2025.svg")

svg(filename = chi2_svg_path, width = 11, height = 12)
print(final_plot)
dev.off()

message("✅ Zapisano plik: ", chi2_svg_path)


# ============================================================
# 💾 Zapis pliku SVG dla OR heatmap
# ============================================================
combined_plot <- wrap_plots(
  overlap_raw_signatures$plot_or,
  overlap_minus_clusters$plot_or,
  overlap_minus_global6$plot_or,
  overlap_minus_global5$plot_or,
  ncol = 2,
  guides = "collect"
) & theme(legend.position = "bottom")

# ============================================================
# 🟩🟥 3️⃣ Utworzenie niestandardowej legendy (puste prostokąty)
# ============================================================

custom_legend <- create_customRect_patch_legend(
  colors = c("#97C426", "#2F4603"),
  labels = c("p < 0.01", "p < 0.0001"),
  spacing = 2,
  box_linewidth = 4
)

# ============================================================
# 🧩 4️⃣ Połączenie wszystkiego (wykresy + legenda)
# ============================================================

final_plot <- wrap_plots(
  combined_plot,
  custom_legend,
  ncol = 1,
  heights = c(10, 1)
)

# ============================================================
# 💾 5️⃣ Podgląd lub zapis
# ============================================================

final_plot

or_svg_path <- file.path(output_dir, "combined_overlap_OR_finalPlot_22.10.2025.svg")

svg(filename = or_svg_path, width = 11, height = 12)
print(final_plot)
dev.off()

message("✅ Zapisano plik: ", or_svg_path)


# ============================================================
# 📁 Ścieżka zapisu wyników overlap analysis
# ============================================================
output_dir <- "/home/mateusz/projects/ifpan-michkor-gr/results_v2/positive_validation_grSignatures_21.10.2025/overlap_chi2"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ============================================================
# 💾 Zapis obiektów overlap (pomiędzy tissue GR signatures)
# ============================================================

saveRDS(
  overlap_raw_signatures,
  file = file.path(
    output_dir,
    "overlap_rawSignatures_tissuesGrSignatures_22.10.2025.rds"
  )
)

saveRDS(
  overlap_minus_clusters,
  file = file.path(
    output_dir,
    "overlap_minusClusters_tissuesGrSignatures_22.10.2025.rds"
  )
)

saveRDS(
  overlap_minus_global6,
  file = file.path(
    output_dir,
    "overlap_minusGlobal6_tissuesGrSignatures_22.10.2025.rds"
  )
)

saveRDS(
  overlap_minus_global5,
  file = file.path(
    output_dir,
    "overlap_minusGlobal5_tissuesGrSignatures_22.10.2025.rds"
  )
)

message("✅ Zapisano wszystkie pliki .rds z wynikami overlap pomiędzy tissue GR signatures:")
message(" - overlap_rawSignatures_tissuesGrSignatures_22.10.2025.rds")
message(" - overlap_minusClusters_tissuesGrSignatures_22.10.2025.rds")
message(" - overlap_minusGlobal6_tissuesGrSignatures_22.10.2025.rds")
message(" - overlap_minusGlobal5_tissuesGrSignatures_22.10.2025.rds")

