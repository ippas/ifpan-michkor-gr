# ============================================================
# 📚 Pakiety
# ============================================================

library(dplyr)
library(purrr)
library(ggplot2)
library(ggsignif)
library(patchwork)
library(openxlsx)
library(stringr)

# ============================================================
# 📁 Ścieżka zapisu wyników
# ============================================================

output_dir <- "/home/mateusz/projects/ifpan-michkor-gr/results_v2/positive_validation_grSignatures_31.10.2025/overlap_chi2/overlap_tissuesGlobalGrSignatures"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ============================================================
# 1️⃣ Raw GR-dependent signatures + global signatures
# ============================================================

label_mapper_raw <- c(
  "NeuralCellsUp" = "NeuralCellsUp",
  "BloodCellsUp" = "BloodCellsUp",
  "LungCellsUp" = "LungCellsUp",
  "NeuralCellsDown" = "NeuralCellsDown",
  "BloodCellsDown" = "BloodCellsDown",
  "LungCellsDown" = "LungCellsDown",
  "global_GR_genes_globalUp5TissuesDerivedCells" = "GlobalUp",
  "global_GR_genes_globalDown5TissuesDerivedCells" = "GlobalDown"
)

overlap_raw_tissues_global <- run_full_overlap_analysis(
  gene_lists = flat_allGrSignatures_31.10.2025[c(
    "NeuralCellsUp", "BloodCellsUp", "LungCellsUp",
    "NeuralCellsDown", "BloodCellsDown", "LungCellsDown",
    "global_GR_genes_globalUp5TissuesDerivedCells",
    "global_GR_genes_globalDown5TissuesDerivedCells"
  )],
  total_genes = hgnc_symbols_vector_v110,
  plot_title_or = "Raw GR-dependent + global signatures — log2(OR)",
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
# 2️⃣ Minus global UP+DOWN (5-tissue criterion) + global signatures
# ============================================================

label_mapper_minusGlobalUpDown <- c(
  "minusGlobalUpDown5TissuesDerivedCells_NeuralCellsUp"   = "NeuralCellsUp",
  "minusGlobalUpDown5TissuesDerivedCells_BloodCellsUp"    = "BloodCellsUp",
  "minusGlobalUpDown5TissuesDerivedCells_LungCellsUp"     = "LungCellsUp",
  "minusGlobalUpDown5TissuesDerivedCells_NeuralCellsDown" = "NeuralCellsDown",
  "minusGlobalUpDown5TissuesDerivedCells_BloodCellsDown"  = "BloodCellsDown",
  "minusGlobalUpDown5TissuesDerivedCells_LungCellsDown"   = "LungCellsDown",
  "global_GR_genes_globalUp5TissuesDerivedCells"          = "GlobalUp",
  "global_GR_genes_globalDown5TissuesDerivedCells"        = "GlobalDown"
)

overlap_minusGlobalUpDown_tissues_global <- run_full_overlap_analysis(
  gene_lists = flat_allGrSignatures_31.10.2025[c(
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
  plot_title_or = "Minus global (UP+DOWN, 5-tissue) + global signatures — log2(OR)",
  plot_title_chi2 = "Minus global (UP+DOWN, 5-tissue) + global signatures — log2(χ² + 1)",
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
  row_labels_map = label_mapper_minusGlobalUpDown,
  col_labels_map = label_mapper_minusGlobalUpDown
)

# ============================================================
# 🧩 Połączenie heatmap — χ²
# ============================================================

combined_plot_chi2 <- wrap_plots(
  overlap_raw_tissues_global$plot_chi2,
  overlap_minusGlobalUpDown_tissues_global$plot_chi2,
  ncol = 2, guides = "collect"
) & theme(legend.position = "bottom")

custom_legend <- create_customRect_patch_legend(
  colors = c("#97C426", "#2F4603"),
  labels = c("p < 0.01", "p < 0.0001"),
  spacing = 4,
  box_linewidth = 4
)

final_plot_chi2 <- wrap_plots(
  combined_plot_chi2, custom_legend, ncol = 1, heights = c(10, 1)
)

# ============================================================
# 💾 Zapis SVG χ²
# ============================================================

chi2_svg_path <- file.path(output_dir, "combined_tissuesGlobalGrSignatures_overlap_chi2_31.10.2025.svg")
svg(filename = chi2_svg_path, width = 14, height = 8)
print(final_plot_chi2)
dev.off()
message("✅ Zapisano plik: ", chi2_svg_path)

# ============================================================
# 🧩 Połączenie heatmap — log2(OR)
# ============================================================

combined_plot_or <- wrap_plots(
  overlap_raw_tissues_global$plot_or,
  overlap_minusGlobalUpDown_tissues_global$plot_or,
  ncol = 2, guides = "collect"
) & theme(legend.position = "bottom")

final_plot_or <- wrap_plots(
  combined_plot_or, custom_legend, ncol = 1, heights = c(10, 1)
)

# ============================================================
# 💾 Zapis SVG OR
# ============================================================

or_svg_path <- file.path(output_dir, "combined_tissuesGlobalGrSignatures_overlap_OR_31.10.2025.svg")
svg(filename = or_svg_path, width = 14, height = 8)
print(final_plot_or)
dev.off()
message("✅ Zapisano plik: ", or_svg_path)

# ============================================================
# 💾 Zapis wyników RDS
# ============================================================

saveRDS(overlap_raw_tissues_global,
        file = file.path(output_dir, "overlap_raw_tissuesGlobalGrSignatures_31.10.2025.rds"))
saveRDS(overlap_minusGlobalUpDown_tissues_global,
        file = file.path(output_dir, "overlap_minusGlobalUpDown_tissuesGlobalGrSignatures_31.10.2025.rds"))

message("✅ Zapisano pliki .rds:")
message(" - overlap_raw_tissuesGlobalGrSignatures_31.10.2025.rds")
message(" - overlap_minusGlobalUpDown_tissuesGlobalGrSignatures_31.10.2025.rds")

# ============================================================
# 1️⃣ Przetwarzanie danych
# ============================================================

df_processed <- overlap_minusGlobalUpDown_tissues_global$processed$original_data$df %>%
  remove_duplicate_pairs() %>%
  mutate(fdr = p.adjust(p_value, method = "BH")) %>%
  select(-any_of("fdr_value")) %>%
  relocate(overlap_genes, .after = last_col()) %>%
  rename(signatureA = Var1, signatureB = Var2) %>%
  mutate(
    across(c(signatureA, signatureB), ~ .x %>%
             str_remove("^minusGlobalUpDown5TissuesDerivedCells_") %>%
             str_replace("^global_GR_genes_globalUp5TissuesDerivedCells$", "globalUp") %>%
             str_replace("^global_GR_genes_globalDown5TissuesDerivedCells$", "globalDown"))
  )

# ============================================================
# 2️⃣ Tworzenie workbooka i arkusza
# ============================================================

wb <- createWorkbook()
addWorksheet(wb, "overlap_results")

# ============================================================
# 3️⃣ Style
# ============================================================

style_num <- createStyle(numFmt = "0.000")
style_sci <- createStyle(numFmt = "0.00E+00")
header_style <- createStyle(textDecoration = "bold", halign = "center")

# ============================================================
# 4️⃣ Zapis danych i formatowanie
# ============================================================

writeData(wb, "overlap_results", df_processed, withFilter = TRUE, headerStyle = header_style)

addStyle(wb, "overlap_results", style = style_sci,
         cols = which(names(df_processed) %in% c("p_value", "fdr")),
         rows = 2:(nrow(df_processed) + 1), gridExpand = TRUE)

addStyle(wb, "overlap_results", style = style_num,
         cols = which(names(df_processed) %in% c("chi2", "gene_overlap_count", "odds_ratio", "log2_odds_ratio")),
         rows = 2:(nrow(df_processed) + 1), gridExpand = TRUE)

# ============================================================
# 5️⃣ Dodatkowe ustawienia
# ============================================================

freezePane(wb, "overlap_results", firstRow = TRUE)
setColWidths(wb, "overlap_results", cols = 1:ncol(df_processed), widths = "auto")

# ============================================================
# 6️⃣ Zapis pliku XLSX
# ============================================================

output_path <- file.path(output_dir, "overlap_grSignaturesMinusGlobalUpDown_5TissuesDerived_31.10.2025.xlsx")
saveWorkbook(wb, output_path, overwrite = TRUE)
message("✅ Zapisano plik: ", output_path)




raw_df <- overlap_minusGlobalUpDown_tissues_global$processed$original_data$df
cat("▶ Surowy wynik: ", nrow(raw_df), " wierszy, ", ncol(raw_df), " kolumn\n")

# ============================================================
# 🧮 Przygotowanie tabeli — minimalne czyszczenie
# ============================================================

df_processed <- raw_df %>%
  rename(signatureA = Var1, signatureB = Var2) %>%
  mutate(
    across(c(signatureA, signatureB), ~ .x %>%
             str_remove("^minusGlobalUpDown5TissuesDerivedCells_") %>%
             str_replace("^global_GR_genes_globalUp5TissuesDerivedCells$", "globalUp") %>%
             str_replace("^global_GR_genes_globalDown5TissuesDerivedCells$", "globalDown"))
  ) %>% remove_duplicate_pairs(col_a = "signatureA", col_b = "signatureB") %>% 
  select(-sig_1, -._key) %>% 
  select(-c(fdr_value, fdr)) %>% 
  mutate(fdr = p.adjust(p_value, method = "BH"))

cat("▶ Po przetworzeniu: ", nrow(df_processed), " wierszy, ", ncol(df_processed), " kolumn\n")

# ============================================================
# 💾 Zapis do pliku Excel z formatowaniem
# ============================================================

output_dir <- "/home/mateusz/projects/ifpan-michkor-gr/results_v2/positive_validation_grSignatures_31.10.2025/overlap_chi2/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

output_file <- file.path(output_dir, "overlap_raw_minusGlobalUpDown_andGlobalSignatures_31.10.2025.xlsx")

# 📘 Utworzenie workbooka
wb <- createWorkbook()
addWorksheet(wb, "raw_overlap_results")

# ✍️ Zapis danych
writeData(wb, "raw_overlap_results", df_processed, withFilter = TRUE)
freezePane(wb, "raw_overlap_results", firstRow = TRUE)

# ============================================================
# 🎨 Formatowanie kolumn
# ============================================================

# 🔹 Styl naukowy (3 miejsca po przecinku)
style_sci <- createStyle(numFmt = "0.000E+00")

# 🔹 Styl zwykły (3 miejsca po przecinku)
style_dec <- createStyle(numFmt = "0.000")

# 🧠 Kolumny do stylu naukowego
sci_cols <- which(names(df_processed) %in% c("p_value", "chi2", "odds_ratio", "fdr"))

# 🧠 Kolumna do stylu zwykłego
dec_cols <- which(names(df_processed) %in% c("log2_odds_ratio"))

# 🎨 Zastosowanie stylów
if (length(sci_cols) > 0) {
  addStyle(
    wb, "raw_overlap_results", style = style_sci,
    rows = 2:(nrow(df_processed) + 1), cols = sci_cols,
    gridExpand = TRUE
  )
}

if (length(dec_cols) > 0) {
  addStyle(
    wb, "raw_overlap_results", style = style_dec,
    rows = 2:(nrow(df_processed) + 1), cols = dec_cols,
    gridExpand = TRUE
  )
}

# 📏 Dostosowanie szerokości kolumn
setColWidths(wb, "raw_overlap_results", cols = 1:ncol(df_processed), widths = "auto")

# 💾 Zapisz plik
saveWorkbook(wb, output_file, overwrite = TRUE)

message("✅ Zapisano surowy wynik do pliku: ", output_file)


message("✅ Zapisano surowy wynik do pliku: ", output_file)

