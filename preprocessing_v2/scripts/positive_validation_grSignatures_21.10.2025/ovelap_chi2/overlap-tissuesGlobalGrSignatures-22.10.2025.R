# ============================================================
# 📚 Pakiety
# ============================================================
library(dplyr)
library(purrr)
library(ggplot2)
library(ggsignif)
library(patchwork)

# ============================================================
# 📁 Ścieżka zapisu wyników
# ============================================================
output_dir <- "/home/mateusz/projects/ifpan-michkor-gr/results_v2/positive_validation_grSignatures_21.10.2025/overlap_chi2/ovelap_tissuesGlobalGrSignatures"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ============================================================
# 1️⃣ Raw GR-dependent signatures + global signatures
# ============================================================
label_mapper_raw <- c(
  "NeuralCellsUp"  = "NeuralCellsUp",
  "BloodCellsUp"   = "BloodCellsUp",
  "LungCellsUp"    = "LungCellsUp",
  "NeuralCellsDown" = "NeuralCellsDown",
  "BloodCellsDown"  = "BloodCellsDown",
  "LungCellsDown"   = "LungCellsDown",
  "global_GR_genes_globalUp5TissuesDerivedCells"   = "GlobalUp",
  "global_GR_genes_globalDown5TissuesDerivedCells" = "GlobalDown"
)

overlap_raw_tissues_global <- run_full_overlap_analysis(
  gene_lists = flat_allGrSignatures_17.10.2025[c(
    "NeuralCellsUp", "BloodCellsUp", "LungCellsUp",
    "NeuralCellsDown", "BloodCellsDown", "LungCellsDown",
    "global_GR_genes_globalUp5TissuesDerivedCells",
    "global_GR_genes_globalDown5TissuesDerivedCells"
  )],
  total_genes = hgnc_symbols_vector_v110,
  
  plot_title_or  = "Raw GR-dependent + global signatures — log2(OR)",
  plot_title_chi2 = "Raw GR-dependent + global signatures — log2(χ² + 1)",
  
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
  p_thresholds_chi2 = c(0.01, 0.0001),
  color_rects_chi2 = c("#97C426", "#2F4603"),
  
  triangle_mode = "upper",
  row_labels_map = label_mapper_raw,
  col_labels_map = label_mapper_raw
)

# ============================================================
# 2️⃣ Minus global (5-tissue criterion) + global signatures
# ============================================================
label_mapper_minusGlobal <- c(
  "minusGlobalUp5TissuesDerivedCells_NeuralCellsUp"     = "NeuralCellsUp",
  "minusGlobalUp5TissuesDerivedCells_BloodCellsUp"      = "BloodCellsUp",
  "minusGlobalUp5TissuesDerivedCells_LungCellsUp"       = "LungCellsUp",
  "minusGlobalDown5TissuesDerivedCells_NeuralCellsDown" = "NeuralCellsDown",
  "minusGlobalDown5TissuesDerivedCells_BloodCellsDown"  = "BloodCellsDown",
  "minusGlobalDown5TissuesDerivedCells_LungCellsDown"   = "LungCellsDown",
  "global_GR_genes_globalUp5TissuesDerivedCells"   = "GlobalUp",
  "global_GR_genes_globalDown5TissuesDerivedCells" = "GlobalDown"
)

overlap_minusGlobal_tissues_global <- run_full_overlap_analysis(
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
  
  plot_title_or  = "Minus global (5-tissue) + global signatures — log2(OR)",
  plot_title_chi2 = "Minus global (5-tissue) + global signatures — log2(χ² + 1)",
  
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  text_size_axis = 14,
  show_dendrograms = FALSE,
  
  # palette_or = c("#07243e", "white", "darkred"), -> darkblue
  palette_or = c("#c6d3e3", "white", "darkred"), # lightblue
  color_scale_range_or = c(-5, 5),
  text_contrast_range_or = c(-30, 4),
  
  palette_chi2 = c("white", "darkred"),
  color_scale_range_chi2 = c(0, 8),
  text_contrast_range_chi2 = c(0, 6),
  p_thresholds_chi2 = c(0.01, 0.0001),
  color_rects_chi2 = c("#97C426", "#2F4603"),
  
  triangle_mode = "upper",
  row_labels_map = label_mapper_minusGlobal,
  col_labels_map = label_mapper_minusGlobal
)

# ============================================================
# 🧩 Połączenie heatmap — χ²
# ============================================================
combined_plot_chi2 <- wrap_plots(
  overlap_raw_tissues_global$plot_chi2,
  overlap_minusGlobal_tissues_global$plot_chi2,
  ncol = 2,
  guides = "collect"
) & theme(legend.position = "bottom")

custom_legend <- create_customRect_patch_legend(
  colors = c("#97C426", "#2F4603"),
  labels = c("p < 0.01", "p < 0.0001"),
  spacing = 4,
  box_linewidth = 4
)

final_plot_chi2 <- wrap_plots(
  combined_plot_chi2,
  custom_legend,
  ncol = 1,
  heights = c(10, 1)
)

# ============================================================
# 💾 Zapis SVG χ²
# ============================================================
chi2_svg_path <- file.path(output_dir, "combined_tissuesGlobalGrSignatures_overlap_chi2_22.10.2025.svg")
svg(filename = chi2_svg_path, width = 14, height = 8)
print(final_plot_chi2)
dev.off()
message("✅ Zapisano plik: ", chi2_svg_path)

# ============================================================
# 🧩 Połączenie heatmap — log2(OR)
# ============================================================
combined_plot_or <- wrap_plots(
  overlap_raw_tissues_global$plot_or,
  overlap_minusGlobal_tissues_global$plot_or,
  ncol = 2,
  guides = "collect"
) & theme(legend.position = "bottom")

final_plot_or <- wrap_plots(
  combined_plot_or,
  custom_legend,
  ncol = 1,
  heights = c(10, 1)
)

# ============================================================
# 💾 Zapis SVG OR
# ============================================================
or_svg_path <- file.path(output_dir, "combined_tissuesGlobalGrSignatures_overlap_OR_22.10.2025.svg")
svg(filename = or_svg_path, width = 14, height = 8)
print(final_plot_or)
dev.off()
message("✅ Zapisano plik: ", or_svg_path)

# ============================================================
# 💾 Zapis wyników RDS
# ============================================================
saveRDS(overlap_raw_tissues_global,
        file = file.path(output_dir, "overlap_raw_tissuesGlobalGrSignatures_22.10.2025.rds"))

saveRDS(overlap_minusGlobal_tissues_global,
        file = file.path(output_dir, "overlap_minusGlobal_tissuesGlobalGrSignatures_22.10.2025.rds"))

message("✅ Zapisano pliki .rds:")
message(" - overlap_raw_tissuesGlobalGrSignatures_22.10.2025.rds")
message(" - overlap_minusGlobal_tissuesGlobalGrSignatures_22.10.2025.rds")

# ============================================================
# 1️⃣ Przetwarzanie danych
# ============================================================
df_processed <- overlap_minusGlobal_tissues_global$processed$original_data$df %>%
  remove_duplicate_pairs() %>%
  mutate(fdr = p.adjust(p_value, method = "BH")) %>%
  select(-any_of("fdr_value")) %>%
  relocate(overlap_genes, .after = last_col()) %>%
  rename(signatureA = Var1, signatureB = Var2) %>%
  mutate(
    across(
      c(signatureA, signatureB),
      ~ .x %>%
        str_remove("^minusGlobal(?:Up|Down)5TissuesDerivedCells_") %>%
        str_replace("^global_GR_genes_globalUp5TissuesDerivedCells$", "globalUp") %>%
        str_replace("^global_GR_genes_globalDown5TissuesDerivedCells$", "globalDown")
    )
  )

# ============================================================
# 2️⃣ Tworzenie workbooka i arkusza
# ============================================================
wb <- createWorkbook()
addWorksheet(wb, "overlap_results")

# ============================================================
# 3️⃣ Style
# ============================================================
style_num <- createStyle(numFmt = "0.000")       # 3 miejsca po przecinku
style_sci <- createStyle(numFmt = "0.00E+00")    # zapis naukowy
header_style <- createStyle(textDecoration = "bold", halign = "center")

# ============================================================
# 4️⃣ Zapis danych i formatowanie
# ============================================================
writeData(
  wb,
  sheet = "overlap_results",
  x = df_processed,
  withFilter = TRUE,
  headerStyle = header_style
)

# format naukowy dla p_value i fdr
addStyle(
  wb, "overlap_results",
  style = style_sci,
  cols = which(names(df_processed) %in% c("p_value", "fdr")),
  rows = 2:(nrow(df_processed) + 1),
  gridExpand = TRUE
)

# 3 miejsca po przecinku dla pozostałych kolumn liczbowych
addStyle(
  wb, "overlap_results",
  style = style_num,
  cols = which(names(df_processed) %in% c("chi2", "gene_overlap_count", "odds_ratio", "log2_odds_ratio")),
  rows = 2:(nrow(df_processed) + 1),
  gridExpand = TRUE
)

# ============================================================
# 5️⃣ Dodatkowe ustawienia
# ============================================================
freezePane(wb, sheet = "overlap_results", firstRow = TRUE)  # blokada pierwszego wiersza
setColWidths(wb, "overlap_results", cols = 1:ncol(df_processed), widths = "auto")  # dopasuj szerokość kolumn

# ============================================================
# 6️⃣ Zapis pliku
# ============================================================
output_path <- "/home/mateusz/projects/ifpan-michkor-gr/results_v2/positive_validation_grSignatures_21.10.2025/overlap_chi2/overlap_grSignaturesMinusGlobal_5TissuesDerived_23.10.2025.xlsx"
saveWorkbook(wb, output_path, overwrite = TRUE)

message("✅ Zapisano plik: ", output_path)

