

# 📋 Nazwy sygnatur do analizy
signatures_to_run <- c(
  
  # ↑ Up-regulated minus global
  "minusGlobalUp5TissuesDerivedCells_BloodCellsUp",
  "minusGlobalUp5TissuesDerivedCells_NeuralCellsUp",
  "minusGlobalUp5TissuesDerivedCells_LungCellsUp",
  
  # ↓ Down-regulated minus global
  "minusGlobalDown5TissuesDerivedCells_BloodCellsDown",
  "minusGlobalDown5TissuesDerivedCells_NeuralCellsDown",
  "minusGlobalDown5TissuesDerivedCells_LungCellsDown",
  
  # 🌍 Global signatures
  "global_GR_genes_globalUp5TissuesDerivedCells",
  "global_GR_genes_globalDown5TissuesDerivedCells"
)

# 🚀 Uruchom Enrichr (ChEA_2022) i przefiltruj wyniki
enrichr_chea_results <- lapply(
  signatures_to_run,
  function(sig_name) {
    cat("▶️ Running Enrichr for:", sig_name, "\n")
    
    result <- run_enrichr(
      gene_list = flat_allGrSignatures_17.10.2025[[sig_name]],
      database = "ChEA_2022"
    )
    
    # 📊 Filtrujemy wyniki
    result %>%
      filter(n_genes > 2) %>%
      filter(Adjusted.P.value < 0.05) %>% 
      select(-c(Old.P.value, Old.Adjusted.P.value))
  }
)

# 📦 Nadaj nazwy elementom listy
names(enrichr_chea_results) <- signatures_to_run


enrichr_chea_results[c(   "minusGlobalDown5TissuesDerivedCells_BloodCellsDown",
                          "minusGlobalDown5TissuesDerivedCells_NeuralCellsDown",
                          "minusGlobalDown5TissuesDerivedCells_LungCellsDown",
                         "global_GR_genes_globalDown5TissuesDerivedCells")] %>% 
  lapply(., function(x){
    x  %>% select(-Genes) %>% filter(grepl("RELB|SUZ12|IRF8|SMRT|THRA|NFE2L2|NRF2", Term))
  })

enrichr_chea_results[c(   "minusGlobalDown5TissuesDerivedCells_BloodCellsDown",
                          "minusGlobalDown5TissuesDerivedCells_NeuralCellsDown",
                          "minusGlobalDown5TissuesDerivedCells_LungCellsDown",
                          "global_GR_genes_globalDown5TissuesDerivedCells")] %>% 
  lapply(., function(x){
    x  %>% select(-Genes) %>% filter(grepl("RELB|RELA|NFKB1|CEBPD|EP300", Term))
  })



enrichr_chea_results$minusGlobalDown5TissuesDerivedCells_BloodCellsDown %>% select

# ============================================================
# 📋 Ustal kolejność i nazwy arkuszy
# ============================================================
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
output_dir <- "/home/mateusz/projects/ifpan-michkor-gr/results_v2/positive_validation_grSignatures_21.10.2025/enrichr_chea/chea_23.10.2025"

# ============================================================
# 🧾 Nazwa pliku wynikowego
# ============================================================
output_file <- file.path(
  output_dir,
  "ChEA2022_grSignaturesGlobal5TissuesDerivedTissuesMinusGlobal_fdr0.05overlap2-23.10.2025.xlsx"
)

# ============================================================
# 🧩 Utwórz workbook
# ============================================================
wb <- createWorkbook()

# ============================================================
# 🎯 Kolumny do formatowania
# ============================================================
numeric_cols <- c("P.value", "Adjusted.P.value", "Odds.Ratio", "Combined.Score")

# ============================================================
# 🧠 Pętla po wszystkich arkuszach
# ============================================================


# 📦 1. Przygotowanie przefiltrowanej listy
filtered_list <- list(
  "NeuralCellsUp"   = enrichr_cellmarker_minusGlobal5TissueDerived$NeuralCells$Up,
  "NeuralCellsDown" = enrichr_cellmarker_minusGlobal5TissueDerived$NeuralCells$Down,
  "BloodCellsUp"    = enrichr_cellmarker_minusGlobal5TissueDerived$BloodCells$Up,
  "BloodCellsDown"  = enrichr_cellmarker_minusGlobal5TissueDerived$BloodCells$Down,
  "LungCellsUp"     = enrichr_cellmarker_minusGlobal5TissueDerived$LungCells$Up,
  "LungCellsDown"   = enrichr_cellmarker_minusGlobal5TissueDerived$LungCells$Down
) %>%
  lapply(function(x) {
    x %>%
      filter(n_genes > 2, Adjusted.P.value < 0.05) %>%
      select(-c(Old.P.value, Old.Adjusted.P.value))
  })

# 📁 2. Ścieżka zapisu
output_dir <- "/home/mateusz/projects/ifpan-michkor-gr/results_v2/positive_validation_grSignatures_21.10.2025/enrichr_cellmaker/cellMarker_23.10.2025"
output_file <- file.path(
  output_dir,
  "cellMarker2024_grSignaturesGlobal5TissuesDerivedTissuesMinusGlobal_fdr0.05overlap3-23.10.2025.xlsx"
)

# 🎨 3. Style do formatowania
style_sci_3dec <- createStyle(numFmt = "0.000E+00")  # notacja naukowa
style_plain_3dec <- createStyle(numFmt = "0.000")     # zwykłe 3 miejsca

# 📘 4. Stwórz workbook
wb <- createWorkbook()

# ✅ Kolejność arkuszy
sheet_order <- c(
  "NeuralCellsUp", "NeuralCellsDown",
  "BloodCellsUp", "BloodCellsDown",
  "LungCellsUp", "LungCellsDown"
)

# 🔁 5. Zapis arkuszy z odpowiednimi stylami
for (sheet in sheet_order) {
  df <- filtered_list[[sheet]]
  addWorksheet(wb, sheet)
  writeData(wb, sheet, df, withFilter = TRUE)
  freezePane(wb, sheet, firstRow = TRUE)
  
  # 🏷️ Dopasuj szerokość kolumny "Term" do zawartości
  if ("Term" %in% colnames(df)) {
    term_col <- which(colnames(df) == "Term")
    setColWidths(wb, sheet, cols = term_col, widths = "auto")
  }
  
  # 🧬 Poszerz kolumnę "Genes"
  if ("Genes" %in% colnames(df)) {
    genes_col <- which(colnames(df) == "Genes")
    setColWidths(wb, sheet, cols = genes_col, widths = "auto")
  }
  
  # 🔬 Styl naukowy (3 miejsca) → P.value, Adjusted.P.value
  sci_cols <- intersect(c("P.value", "Adjusted.P.value"), colnames(df))
  for (col_name in sci_cols) {
    col_index <- which(colnames(df) == col_name)
    addStyle(wb, sheet, style = style_sci_3dec,
             rows = 2:(nrow(df) + 1), cols = col_index,
             gridExpand = TRUE)
  }
  
  # 📊 Styl zwykły (3 miejsca) → Odds.Ratio, Combined.Score
  plain_cols <- intersect(c("Odds.Ratio", "Combined.Score"), colnames(df))
  for (col_name in plain_cols) {
    col_index <- which(colnames(df) == col_name)
    addStyle(wb, sheet, style = style_plain_3dec,
             rows = 2:(nrow(df) + 1), cols = col_index,
             gridExpand = TRUE)
  }
}

# 💾 6. Zapisz workbook
saveWorkbook(wb, output_file, overwrite = TRUE)

message("✅ Plik zapisany: ", output_file)
