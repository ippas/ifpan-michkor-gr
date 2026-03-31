# ============================================================
# 📦 Wczytanie pakietów
# ============================================================
library(dplyr)
library(purrr)
library(openxlsx)

# ============================================================
# 📋 Nazwy sygnatur do analizy (nowe listy: minusGlobalUpDown5TissuesDerivedCells)
# ============================================================
signatures_to_run <- c(
  # ↑ Up-regulated
  "minusGlobalUpDown5TissuesDerivedCells_BloodCellsUp",
  "minusGlobalUpDown5TissuesDerivedCells_NeuralCellsUp",
  "minusGlobalUpDown5TissuesDerivedCells_LungCellsUp",
  
  # ↓ Down-regulated
  "minusGlobalUpDown5TissuesDerivedCells_BloodCellsDown",
  "minusGlobalUpDown5TissuesDerivedCells_NeuralCellsDown",
  "minusGlobalUpDown5TissuesDerivedCells_LungCellsDown",
  
  # 🌍 Global signatures
  "global_GR_genes_globalUp5TissuesDerivedCells",
  "global_GR_genes_globalDown5TissuesDerivedCells"
)

# ============================================================
# 🚀 Uruchom Enrichr (ChEA_2022)
# ============================================================
enrichr_chea_results <- lapply(
  signatures_to_run,
  function(sig_name) {
    cat("▶️ Running Enrichr for:", sig_name, "\n")
    
    result <- run_enrichr(
      gene_list = flat_allGrSignatures_31.10.2025[[sig_name]],
      database = "ChEA_2022"
    )
    
    result %>%
      filter(n_genes > 2) %>%
      filter(Adjusted.P.value < 0.05) %>%
      select(-c(Old.P.value, Old.Adjusted.P.value))
  }
)

# 📦 Nadaj nazwy elementom listy
names(enrichr_chea_results) <- signatures_to_run

# ============================================================
# 🔍 Przykładowe filtrowanie wyników — RELB, SUZ12, IRF8, itp.
# ============================================================
enrichr_chea_results[c(
  "minusGlobalUpDown5TissuesDerivedCells_BloodCellsDown",
  "minusGlobalUpDown5TissuesDerivedCells_NeuralCellsDown",
  "minusGlobalUpDown5TissuesDerivedCells_LungCellsDown",
  "global_GR_genes_globalDown5TissuesDerivedCells"
)] %>%
  lapply(., function(x) {
    x %>% select(-Genes) %>%
      filter(grepl("RELB|SUZ12|IRF8|SMRT|THRA|NFE2L2|NRF2", Term))
  })

enrichr_chea_results[c(
  "minusGlobalUpDown5TissuesDerivedCells_BloodCellsDown",
  "minusGlobalUpDown5TissuesDerivedCells_NeuralCellsDown",
  "minusGlobalUpDown5TissuesDerivedCells_LungCellsDown",
  "global_GR_genes_globalDown5TissuesDerivedCells"
)] %>%
  lapply(., function(x) {
    x %>% select(-Genes) %>%
      filter(grepl("RELB|RELA|NFKB1|CEBPD|EP300", Term))
  })

# ============================================================
# 📋 Ustal kolejność i nazwy arkuszy
# ============================================================
sheet_order <- c(
  "global_GR_genes_globalUp5TissuesDerivedCells",
  "global_GR_genes_globalDown5TissuesDerivedCells",
  "minusGlobalUpDown5TissuesDerivedCells_NeuralCellsUp",
  "minusGlobalUpDown5TissuesDerivedCells_NeuralCellsDown",
  "minusGlobalUpDown5TissuesDerivedCells_BloodCellsUp",
  "minusGlobalUpDown5TissuesDerivedCells_BloodCellsDown",
  "minusGlobalUpDown5TissuesDerivedCells_LungCellsUp",
  "minusGlobalUpDown5TissuesDerivedCells_LungCellsDown"
)

sheet_names <- c(
  "globalGrUp_5TissuesDerived",
  "globalGrDown_5TissuesDerived",
  "NeuralCellsUp",
  "NeuralCellsDown",
  "BloodCellsUp",
  "BloodCellsDown",
  "LungCellsUp",
  "LungCellsDown"
)

# ============================================================
# 📂 Ścieżka do folderu zapisu
# ============================================================
output_dir <- "/home/mateusz/projects/ifpan-michkor-gr/results_v2/positive_validation_grSignatures_31.10.2025/enrichr_chea/chea_31.10.2025"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 🧾 Nazwa pliku wynikowego
# ============================================================
output_file <- file.path(
  output_dir,
  "ChEA2022_grSignaturesGlobal5TissuesDerivedTissuesMinusGlobalUpDown_fdr0.05overlap2-31.10.2025.xlsx"
)

# ============================================================
# 🎨 Style do formatowania
# ============================================================
style_sci_3dec <- createStyle(numFmt = "0.000E+00")  # notacja naukowa
style_plain_3dec <- createStyle(numFmt = "0.000")     # zwykłe 3 miejsca

# ============================================================
# 🧩 Utwórz workbook
# ============================================================
wb <- createWorkbook()

# ============================================================
# 🧠 Pętla po wszystkich arkuszach
# ============================================================
for (i in seq_along(sheet_order)) {
  sig_name <- sheet_order[i]
  sheet_name <- sheet_names[i]
  df <- enrichr_chea_results[[sig_name]]
  
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, df, withFilter = TRUE)
  freezePane(wb, sheet_name, firstRow = TRUE)
  
  # 🏷️ Dopasuj szerokości kolumn
  if ("Term" %in% colnames(df)) {
    term_col <- which(colnames(df) == "Term")
    setColWidths(wb, sheet_name, cols = term_col, widths = "auto")
  }
  if ("Genes" %in% colnames(df)) {
    genes_col <- which(colnames(df) == "Genes")
    setColWidths(wb, sheet_name, cols = genes_col, widths = "auto")
  }
  
  # 🔬 Styl naukowy dla p-value i adjusted p-value
  sci_cols <- intersect(c("P.value", "Adjusted.P.value"), colnames(df))
  for (col_name in sci_cols) {
    col_index <- which(colnames(df) == col_name)
    addStyle(
      wb, sheet_name, style = style_sci_3dec,
      rows = 2:(nrow(df) + 1), cols = col_index,
      gridExpand = TRUE
    )
  }
  
  # 📊 Styl zwykły (3 miejsca po przecinku)
  plain_cols <- intersect(c("Odds.Ratio", "Combined.Score"), colnames(df))
  for (col_name in plain_cols) {
    col_index <- which(colnames(df) == col_name)
    addStyle(
      wb, sheet_name, style = style_plain_3dec,
      rows = 2:(nrow(df) + 1), cols = col_index,
      gridExpand = TRUE
    )
  }
}

# ============================================================
# 💾 Zapisz workbook
# ============================================================
saveWorkbook(wb, output_file, overwrite = TRUE)

message("✅ Zapisano plik: ", output_file)
