
# 📋 Nazwy arkuszy i kolumn
sheet_order <- c(
  "global_GR_genes_globalUp5TissuesDerivedCells",
  "global_GR_genes_globalDown5TissuesDerivedCells",
  "minusGlobalUp5TissuesDerivedCells_NeuralCellsUp",
  "minusGlobalDown5TissuesDerivedCells_NeuralCellsDown",
  "minusGlobalUp5TissuesDerivedCells_BloodCellsUp",
  "minusGlobalDown5TissuesDerivedCells_BloodCellsDown",
  "minusGlobalUp5TissuesDerivedCells_LungCellsUp",
  "minusGlobalDown5TissuesDerivedCells_LungCellsDown"
)

column_names <- c(
  "globalGrUp_5TissuesDerived",
  "globalGrDown_5TissuesDerived",
  "NeuralCellsUp",
  "NeuralCellsDown",
  "BloodCellsUp",
  "BloodCellsDown",
  "LungCellsUp",
  "LungCellsDown"
)

# 📂 Ścieżka zapisu
output_dir <- "/home/mateusz/projects/ifpan-michkor-gr/results_v2/GR_signatures/GR_signatures_24.10.2025"
output_file <- file.path(output_dir, "grSignaturures_global5TsissuesDerived_tissuesMinusGlobal_24.10.2025.xlsx")

# 🧩 Utworzenie jednej tabeli z różnymi długościami kolumn
gene_lists <- lapply(sheet_order, function(x) flat_allGrSignatures_17.10.2025[[x]])
max_len <- max(sapply(gene_lists, length))
gene_df <- purrr::map_dfc(gene_lists, ~ {
  vec <- .x
  length(vec) <- max_len  # uzupełnij NA
  tibble(vec)
})
colnames(gene_df) <- column_names

# 💾 Zapis do XLSX
wb <- createWorkbook()
addWorksheet(wb, "GR_signatures")

# Zapis danych
writeData(wb, "GR_signatures", gene_df, headerStyle = createStyle(textDecoration = "bold"))

# Zablokuj pierwszy wiersz
freezePane(wb, "GR_signatures", firstActiveRow = 2)

# Dopasuj szerokość kolumn
setColWidths(wb, "GR_signatures", cols = 1:ncol(gene_df), widths = "auto")

# Zapisz plik
saveWorkbook(wb, output_file, overwrite = TRUE)

message("✅ Zapisano plik: ", output_file)