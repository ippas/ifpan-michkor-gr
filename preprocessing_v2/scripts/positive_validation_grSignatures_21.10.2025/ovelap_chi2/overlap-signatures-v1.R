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
  
  plot_title_or = "Raw GR-dependent signatures — log2(OR)",
  plot_title_chi2 = "Raw GR-dependent signatures — log2(χ² + 1)",
  
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  text_size_axis = 14,
  show_dendrograms = FALSE,
  
  palette_or = c("#07243e", "white", "darkred"),
  color_scale_range_or = c(-5, 5),
  text_contrast_range_or = c(-3, 4.9),
  
  palette_chi2 = c("white", "darkred"),
  color_scale_range_chi2 = c(0, 10),
  text_contrast_range_chi2 = c(0, 4),
  p_thresholds_chi2 = c(0.01, 0.0001),
  color_rects_chi2 = c("#97C426", "#2F4603"),
  
  triangle_mode = "upper"
)

# ============================================================
# 2️⃣ Minus clusters (KPOD/BMC removed)
# ============================================================
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
  color_scale_range_chi2 = c(0, 12),
  text_contrast_range_chi2 = c(0, 4),
  p_thresholds_chi2 = c(0.01, 0.0001),
  color_rects_chi2 = c("#97C426", "#2F4603"),
  
  triangle_mode = "upper"
)

# ============================================================
# 3️⃣ Minus global (4 tissues)
# ============================================================
overlap_minus_global4 <- run_full_overlap_analysis(
  gene_lists = flat_allGrSignatures_17.10.2025[c(
    "minusGlobalUp4TissuesDerivedCells_NeuralCellsUp",
    "minusGlobalUp4TissuesDerivedCells_BloodCellsUp",
    "minusGlobalUp4TissuesDerivedCells_LungCellsUp",
    "minusGlobalDown4TissuesDerivedCells_NeuralCellsDown",
    "minusGlobalDown4TissuesDerivedCells_BloodCellsDown",
    "minusGlobalDown4TissuesDerivedCells_LungCellsDown"
  )],
  total_genes = hgnc_symbols_vector_v110,
  
  plot_title_or = "Minus global (4-tissue criterion) — log2(OR)",
  plot_title_chi2 = "Minus global (4-tissue criterion) — log2(χ² + 1)",
  
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  text_size_axis = 14,
  show_dendrograms = FALSE,
  
  palette_or = c("#07243e", "white", "darkred"),
  color_scale_range_or = c(-5, 5),
  text_contrast_range_or = c(-3, 4.9),
  
  palette_chi2 = c("white", "darkred"),
  color_scale_range_chi2 = c(0, 12),
  text_contrast_range_chi2 = c(0, 4),
  p_thresholds_chi2 = c(0.01, 0.0001),
  color_rects_chi2 = c("#97C426", "#2F4603"),
  
  triangle_mode = "upper"
)

# ============================================================
# 4️⃣ Minus global (5 tissues)
# ============================================================
# ============================================================
# Mapper etykiet — usuwa prefiksy "minusGlobalUp/DownX..."
# ============================================================
label_mapper <- c(
  "minusGlobalUp5TissuesDerivedCells_NeuralCellsUp"    = "NeuralCellsUp",
  "minusGlobalUp5TissuesDerivedCells_BloodCellsUp"     = "BloodCellsUp",
  "minusGlobalUp5TissuesDerivedCells_LungCellsUp"      = "LungCellsUp",
  "minusGlobalDown5TissuesDerivedCells_NeuralCellsDown" = "NeuralCellsDown",
  "minusGlobalDown5TissuesDerivedCells_BloodCellsDown"  = "BloodCellsDown",
  "minusGlobalDown5TissuesDerivedCells_LungCellsDown"   = "LungCellsDown"
)


# ============================================================
# Analiza: minus global (5 tissues)
# ============================================================
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
  color_scale_range_chi2 = c(0, 10),
  text_contrast_range_chi2 = c(0, 4),
  p_thresholds_chi2 = c(0.01, 0.0001),
  # p_thresholds
  color_rects_chi2 = c("#97C426", "#2F4603"),
  
  triangle_mode = "upper",
  
  # # 👇 przekazanie mappera do obu heatmap
  row_labels_map = label_mapper,
  col_labels_map = label_mapper
)


  overlap_minus_global5$plot_chi2

