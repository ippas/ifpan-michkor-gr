# ============================================================
# 📚 Pakiety
# ============================================================

library(dplyr)
library(purrr)
library(ggplot2)
library(patchwork)
library(openxlsx)
library(stringr)

# ============================================================
# 📁 Ścieżka zapisu wyników
# ============================================================

output_dir <- "/home/mateusz/projects/ifpan-michkor-gr/results_v2/positive_validation_grSignatures_18.11.2025/overlap_chi2/overlap_tissuesGlobalGrSignatures"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ============================================================
# 🏷 MAPOWANIA LABELI (oś + nazwy globalUp/globalDown)
# ============================================================

label_mapper_raw <- c(
  "NeuralCellsUp"   = "NeuralCellsUp",
  "BloodCellsUp"    = "BloodCellsUp",
  "LungCellsUp"     = "LungCellsUp",
  "NeuralCellsDown" = "NeuralCellsDown",
  "BloodCellsDown"  = "BloodCellsDown",
  "LungCellsDown"   = "LungCellsDown",
  "global_GR_genes_globalUp5TissuesDerivedCells"   = "GlobalUp5Tissues",
  "global_GR_genes_globalDown5TissuesDerivedCells" = "GlobalDown5Tissues"
)

label_mapper_minusGlobalUpDown4 <- c(
  "minusGlobalUpDown4TissuesDerivedCells_NeuralCellsUp"   = "NeuralCellsUp",
  "minusGlobalUpDown4TissuesDerivedCells_BloodCellsUp"    = "BloodCellsUp",
  "minusGlobalUpDown4TissuesDerivedCells_LungCellsUp"     = "LungCellsUp",
  "minusGlobalUpDown4TissuesDerivedCells_NeuralCellsDown" = "NeuralCellsDown",
  "minusGlobalUpDown4TissuesDerivedCells_BloodCellsDown"  = "BloodCellsDown",
  "minusGlobalUpDown4TissuesDerivedCells_LungCellsDown"   = "LungCellsDown",
  "global_GR_genes_globalUp4TissuesDerivedCells"          = "GlobalUp4Tissues",
  "global_GR_genes_globalDown4TissuesDerivedCells"        = "GlobalDown4Tissues"
)

label_mapper_minusGlobalUpDown5 <- c(
  "minusGlobalUpDown5TissuesDerivedCells_NeuralCellsUp"   = "NeuralCellsUp",
  "minusGlobalUpDown5TissuesDerivedCells_BloodCellsUp"    = "BloodCellsUp",
  "minusGlobalUpDown5TissuesDerivedCells_LungCellsUp"     = "LungCellsUp",
  "minusGlobalUpDown5TissuesDerivedCells_NeuralCellsDown" = "NeuralCellsDown",
  "minusGlobalUpDown5TissuesDerivedCells_BloodCellsDown"  = "BloodCellsDown",
  "minusGlobalUpDown5TissuesDerivedCells_LungCellsDown"   = "LungCellsDown",
  "global_GR_genes_globalUp5TissuesDerivedCells"          = "GlobalUp5Tissues",
  "global_GR_genes_globalDown5TissuesDerivedCells"        = "GlobalDown5Tissues"
)

label_mapper_minusGlobalUpDown6 <- c(
  "minusGlobalUpDown6TissuesDerivedCells_NeuralCellsUp"   = "NeuralCellsUp",
  "minusGlobalUpDown6TissuesDerivedCells_BloodCellsUp"    = "BloodCellsUp",
  "minusGlobalUpDown6TissuesDerivedCells_LungCellsUp"     = "LungCellsUp",
  "minusGlobalUpDown6TissuesDerivedCells_NeuralCellsDown" = "NeuralCellsDown",
  "minusGlobalUpDown6TissuesDerivedCells_BloodCellsDown"  = "BloodCellsDown",
  "minusGlobalUpDown6TissuesDerivedCells_LungCellsDown"   = "LungCellsDown",
  "global_GR_genes_globalUp6TissuesDerivedCells"          = "GlobalUp6Tissues",
  "global_GR_genes_globalDown6TissuesDerivedCells"        = "GlobalDown6Tissues"
)

# ============================================================
# 🔥 1) ANALIZY OVERLAP – RAW + 4/5/6 TISSUES
#    (palety takie jak miałeś – NIC nie zmieniam)
# ============================================================

# RAW (tissue signatures + global 5-tissue)
overlap_raw_tissues_global <- run_full_overlap_analysis(
  gene_lists = flat_allGrSignatures_18.11.2025[c(
    "NeuralCellsUp", "BloodCellsUp", "LungCellsUp",
    "NeuralCellsDown", "BloodCellsDown", "LungCellsDown",
    "global_GR_genes_globalUp5TissuesDerivedCells",
    "global_GR_genes_globalDown5TissuesDerivedCells"
  )],
  total_genes = hgnc_symbols_vector_v110,
  plot_title_or   = "RAW GR-dependent + GlobalUp5Tissues / GlobalDown5Tissues — log2(OR)",
  plot_title_chi2 = "RAW GR-dependent + GlobalUp5Tissues / GlobalDown5Tissues — log2(χ² + 1)",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  text_size_axis = 14,
  show_dendrograms = FALSE,
  palette_or = c("#c6d3e3", "white", "darkred"),
  color_scale_range_or = c(-5, 5),
  text_contrast_range_or = c(-30, 4),
  palette_chi2 = c("white", "darkred"),
  color_scale_range_chi2 = c(0, 8),
  text_contrast_range_chi2 = c(0, 6),
  p_thresholds_or = c(0.05, 0.01),
  color_rects_or = c("#97C426", "#2F4603"),
  p_thresholds_chi2 = c(0.05, 0.01),
  color_rects_chi2 = c("#97C426", "#2F4603"),
  triangle_mode = "upper",
  row_labels_map = label_mapper_raw,
  col_labels_map = label_mapper_raw
)

# minusGlobalUpDown – 4 tissues
overlap_minusGlobalUpDown4_tissues_global <- run_full_overlap_analysis(
  gene_lists = flat_allGrSignatures_18.11.2025[c(
    "minusGlobalUpDown4TissuesDerivedCells_NeuralCellsUp",
    "minusGlobalUpDown4TissuesDerivedCells_BloodCellsUp",
    "minusGlobalUpDown4TissuesDerivedCells_LungCellsUp",
    "minusGlobalUpDown4TissuesDerivedCells_NeuralCellsDown",
    "minusGlobalUpDown4TissuesDerivedCells_BloodCellsDown",
    "minusGlobalUpDown4TissuesDerivedCells_LungCellsDown",
    "global_GR_genes_globalUp4TissuesDerivedCells",
    "global_GR_genes_globalDown4TissuesDerivedCells"
  )],
  total_genes = hgnc_symbols_vector_v110,
  plot_title_or   = "MinusGlobalUpDown4Tissues + GlobalUp4Tissues / GlobalDown4Tissues — log2(OR)",
  plot_title_chi2 = "MinusGlobalUpDown4Tissues + GlobalUp4Tissues / GlobalDown4Tissues — log2(χ² + 1)",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  text_size_axis = 14,
  show_dendrograms = FALSE,
  palette_or = c("#c6d3e3", "white", "darkred"),
  color_scale_range_or = c(-5, 5),
  text_contrast_range_or = c(-30, 4),
  palette_chi2 = c("white", "darkred"),
  color_scale_range_chi2 = c(0, 8),
  text_contrast_range_chi2 = c(0, 6),
  p_thresholds_or = c(0.05, 0.01),
  color_rects_or = c("#97C426", "#2F4603"),
  p_thresholds_chi2 = c(0.05, 0.01),
  color_rects_chi2 = c("#97C426", "#2F4603"),
  triangle_mode = "upper",
  row_labels_map = label_mapper_minusGlobalUpDown4,
  col_labels_map = label_mapper_minusGlobalUpDown4
)

# minusGlobalUpDown – 5 tissues
overlap_minusGlobalUpDown5_tissues_global <- run_full_overlap_analysis(
  gene_lists = flat_allGrSignatures_18.11.2025[c(
    "minusGlobalUpDown5TissuesDerivedCells_NeuralCellsUp",
    "minusGlobalUpDown5TissuesDerivedCells_BloodCellsUp",
    "minusGlobalUpDown5TissuesDerivedCells_LungCellsUp",
    "minusGlobalUpDown5TissuesDerivedCells_NeuralCellsDown",
    "minusGlobalUpDown5TissuesDerivedCells_BloodCellsDown",
    "minusGlobalUpDown5TissuesDerivedCells_LungCellsDown",
    "global_GR_genes_globalUp5TissuesDerivedCells",
    "global_GR_genes_globalDown5TissuesDerivedCells"
  )],
  total_genes = hgnc_symbols_vector_v110,
  plot_title_or   = "MinusGlobalUpDown5Tissues + GlobalUp5Tissues / GlobalDown5Tissues — log2(OR)",
  plot_title_chi2 = "MinusGlobalUpDown5Tissues + GlobalUp5Tissues / GlobalDown5Tissues — log2(χ² + 1)",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  text_size_axis = 14,
  show_dendrograms = FALSE,
  palette_or = c("#c6d3e3", "white", "darkred"),
  color_scale_range_or = c(-5, 5),
  text_contrast_range_or = c(-30, 4),
  palette_chi2 = c("white", "darkred"),
  color_scale_range_chi2 = c(0, 8),
  text_contrast_range_chi2 = c(0, 6),
  p_thresholds_or = c(0.05, 0.01),
  color_rects_or = c("#97C426", "#2F4603"),
  p_thresholds_chi2 = c(0.05, 0.01),
  color_rects_chi2 = c("#97C426", "#2F4603"),
  triangle_mode = "upper",
  row_labels_map = label_mapper_minusGlobalUpDown5,
  col_labels_map = label_mapper_minusGlobalUpDown5
)

# minusGlobalUpDown – 6 tissues
overlap_minusGlobalUpDown6_tissues_global <- run_full_overlap_analysis(
  gene_lists = flat_allGrSignatures_18.11.2025[c(
    "minusGlobalUpDown6TissuesDerivedCells_NeuralCellsUp",
    "minusGlobalUpDown6TissuesDerivedCells_BloodCellsUp",
    "minusGlobalUpDown6TissuesDerivedCells_LungCellsUp",
    "minusGlobalUpDown6TissuesDerivedCells_NeuralCellsDown",
    "minusGlobalUpDown6TissuesDerivedCells_BloodCellsDown",
    "minusGlobalUpDown6TissuesDerivedCells_LungCellsDown",
    "global_GR_genes_globalUp6TissuesDerivedCells",
    "global_GR_genes_globalDown6TissuesDerivedCells"
  )],
  total_genes = hgnc_symbols_vector_v110,
  plot_title_or   = "MinusGlobalUpDown6Tissues + GlobalUp6Tissues / GlobalDown6Tissues — log2(OR)",
  plot_title_chi2 = "MinusGlobalUpDown6Tissues + GlobalUp6Tissues / GlobalDown6Tissues — log2(χ² + 1)",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  text_size_axis = 14,
  show_dendrograms = FALSE,
  palette_or = c("#c6d3e3", "white", "darkred"),
  color_scale_range_or = c(-5, 5),
  text_contrast_range_or = c(-30, 4),
  palette_chi2 = c("white", "darkred"),
  color_scale_range_chi2 = c(0, 8),
  text_contrast_range_chi2 = c(0, 6),
  p_thresholds_or = c(0.05, 0.01),
  color_rects_or = c("#97C426", "#2F4603"),
  p_thresholds_chi2 = c(0.05, 0.01),
  color_rects_chi2 = c("#97C426", "#2F4603"),
  triangle_mode = "upper",
  row_labels_map = label_mapper_minusGlobalUpDown6,
  col_labels_map = label_mapper_minusGlobalUpDown6

)
# ============================================================
# 🔗 2) POŁĄCZONE HEATMAPY – χ²
# ============================================================

# tutaj NIE ruszam palet – biorę gotowe plot_chi2 z obiektów
combined_plot_chi2 <- wrap_plots(
  overlap_raw_tissues_global$plot_chi2 +
    ggtitle("RAW GR-dependent + GlobalUp5Tissues / GlobalDown5Tissues"),
  
  overlap_minusGlobalUpDown4_tissues_global$plot_chi2 +
    ggtitle("MinusGlobalUpDown4Tissues + GlobalUp4Tissues / GlobalDown4Tissues"),
  
  overlap_minusGlobalUpDown5_tissues_global$plot_chi2 +
    ggtitle("MinusGlobalUpDown5Tissues + GlobalUp5Tissues / GlobalDown5Tissues"),
  
  overlap_minusGlobalUpDown6_tissues_global$plot_chi2 +
    ggtitle("MinusGlobalUpDown6Tissues + GlobalUp6Tissues / GlobalDown6Tissues"),
  
  ncol = 2,
  guides = "collect"
) & theme(
  legend.position = "bottom",
  plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
)

custom_legend <- create_customRect_patch_legend(
  colors = c("#97C426", "#2F4603"),
  labels = c("p < 0.05", "p < 0.01"),
  spacing = 4,
  box_linewidth = 4
)

final_plot_chi2 <- wrap_plots(
  combined_plot_chi2,
  custom_legend,
  ncol = 1,
  heights = c(10, 1)
)

chi2_svg_path <- file.path(output_dir, "combined_tissuesGlobalGrSignatures_overlap_chi2_18.11.2025.svg")
svg(filename = chi2_svg_path, width = 26, height = 26)
print(final_plot_chi2)
dev.off()
message("✅ Zapisano plik χ² SVG: ", chi2_svg_path)

# ============================================================
# 🔗 3) POŁĄCZONE HEATMAPY – log2(OR)
# ============================================================

combined_plot_or <- wrap_plots(
  overlap_raw_tissues_global$plot_or +
    ggtitle("RAW GR-dependent + GlobalUp5Tissues / GlobalDown5Tissues"),
  
  overlap_minusGlobalUpDown4_tissues_global$plot_or +
    ggtitle("MinusGlobalUpDown4Tissues + GlobalUp4Tissues / GlobalDown4Tissues"),
  
  overlap_minusGlobalUpDown5_tissues_global$plot_or +
    ggtitle("MinusGlobalUpDown5Tissues + GlobalUp5Tissues / GlobalDown5Tissues"),
  
  overlap_minusGlobalUpDown6_tissues_global$plot_or +
    ggtitle("MinusGlobalUpDown6Tissues + GlobalUp6Tissues / GlobalDown6Tissues"),
  
  ncol = 2,
  guides = "collect"
) & theme(
  legend.position = "bottom",
  plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
)

final_plot_or <- wrap_plots(
  combined_plot_or,
  custom_legend,
  ncol = 1,
  heights = c(10, 1)
)

or_svg_path <- file.path(output_dir, "combined_tissuesGlobalGrSignatures_overlap_OR_18.11.2025.svg")
svg(filename = or_svg_path, width = 14, height = 14)
print(final_plot_or)
dev.off()
message("✅ Zapisano plik OR SVG: ", or_svg_path)

# ============================================================
# 💾 4) ZAPIS OBIEKTÓW RDS
# ============================================================

saveRDS(
  overlap_raw_tissues_global,
  file = file.path(output_dir, "overlap_raw_tissuesGlobalGrSignatures_18.11.2025.rds")
)
saveRDS(
  overlap_minusGlobalUpDown4_tissues_global,
  file = file.path(output_dir, "overlap_minusGlobalUpDown4TissuesDerived_tissuesGlobalGrSignatures_18.11.2025.rds")
)
saveRDS(
  overlap_minusGlobalUpDown5_tissues_global,
  file = file.path(output_dir, "overlap_minusGlobalUpDown5TissuesDerived_tissuesGlobalGrSignatures_18.11.2025.rds")
)
saveRDS(
  overlap_minusGlobalUpDown6_tissues_global,
  file = file.path(output_dir, "overlap_minusGlobalUpDown6TissuesDerived_tissuesGlobalGrSignatures_18.11.2025.rds")
)

message("✅ Zapisano obiekty .rds do: ", output_dir)

# ===============================================
# 4 pliki XLSX
# ===============================================
process_and_export <- function(df_object, file_name, output_dir) {
  
  df_clean <- df_object$processed$original_data$df %>%
    remove_duplicate_pairs(col_a = "Var1", col_b = "Var2") %>%
    select(-c(._key, fdr, fdr_value, sig_1)) %>%
    mutate(
      p_value = formatC(p_value, format = "e", digits = 3),
      chi2 = round(chi2, 3),
      odds_ratio = round(odds_ratio, 3),
      log2_odds_ratio = round(log2_odds_ratio, 3)
    )
  
  wb <- createWorkbook()
  addWorksheet(wb, sheetName = "results")
  writeData(wb, sheet = "results", x = df_clean, rowNames = FALSE)
  freezePane(wb, sheet = "results", firstRow = TRUE)
  setColWidths(wb, sheet = "results", cols = 1:ncol(df_clean), widths = "auto")
  
  saveWorkbook(
    wb,
    file = file.path(output_dir, file_name),
    overwrite = TRUE
  )
}



process_and_export(overlap_raw_tissues_global,
                   "overlap_raw_tissuesGlobal_18.11.2025.xlsx",
                   output_dir)

process_and_export(overlap_minusGlobalUpDown4_tissues_global,
                   "overlap_minusGlobalUpDown4Tissues_18.11.2025.xlsx",
                   output_dir)

process_and_export(overlap_minusGlobalUpDown5_tissues_global,
                   "overlap_minusGlobalUpDown5Tissues_18.11.2025.xlsx",
                   output_dir)

process_and_export(overlap_minusGlobalUpDown6_tissues_global,
                   "overlap_minusGlobalUpDown6Tissues_18.11.2025.xlsx",
                   output_dir)

