# ##############################################################################
# ---- summary globalUpDown4TissuesDerived ----
# ##############################################################################
# global
flat_allGrSignatures_18.11.2025$global_GR_genes_globalUp4TissuesDerivedCells %>% length()
flat_allGrSignatures_18.11.2025$global_GR_genes_globalDown4TissuesDerivedCells %>% length()


intersect(flat_allGrSignatures_18.11.2025$global_GR_genes_globalUp4TissuesDerivedCells,
          flat_allGrSignatures_18.11.2025$global_GR_genes_globalDown4TissuesDerivedCells
) %>% length()


# blood
flat_allGrSignatures_18.11.2025$minusGlobalUpDown4TissuesDerivedCells_BloodCellsUp %>% length()
flat_allGrSignatures_18.11.2025$minusGlobalUpDown4TissuesDerivedCells_BloodCellsDown %>% length()


intersect(flat_allGrSignatures_18.11.2025$minusGlobalUpDown4TissuesDerivedCells_BloodCellsUp,
          flat_allGrSignatures_18.11.2025$minusGlobalUpDown4TissuesDerivedCells_BloodCellsDown
          )


# neural
flat_allGrSignatures_18.11.2025$minusGlobalUpDown4TissuesDerivedCells_NeuralCellsUp %>% length()
flat_allGrSignatures_18.11.2025$minusGlobalUpDown4TissuesDerivedCells_NeuralCellsDown %>% length()

intersect(flat_allGrSignatures_18.11.2025$minusGlobalUpDown4TissuesDerivedCells_NeuralCellsUp,
          flat_allGrSignatures_18.11.2025$minusGlobalUpDown4TissuesDerivedCells_NeuralCellsDown
)

# lung
flat_allGrSignatures_18.11.2025$minusGlobalUpDown4TissuesDerivedCells_LungCellsUp %>% length()
flat_allGrSignatures_18.11.2025$minusGlobalUpDown4TissuesDerivedCells_LungCellsDown %>% length()

intersect(flat_allGrSignatures_18.11.2025$minusGlobalUpDown4TissuesDerivedCells_LungCellsUp,
          flat_allGrSignatures_18.11.2025$minusGlobalUpDown4TissuesDerivedCells_LungCellsDown
)

# ##############################################################################
# ---- summary globalUpDown6TissuesDerived ----
# ##############################################################################
# global
flat_allGrSignatures_18.11.2025$global_GR_genes_globalUp6TissuesDerivedCells %>% length()
flat_allGrSignatures_18.11.2025$global_GR_genes_globalDown6TissuesDerivedCells %>% length()

# blood
flat_allGrSignatures_18.11.2025$minusGlobalUpDown6TissuesDerivedCells_BloodCellsUp %>% length()
flat_allGrSignatures_18.11.2025$minusGlobalUpDown6TissuesDerivedCells_BloodCellsDown %>% length()


intersect(flat_allGrSignatures_18.11.2025$minusGlobalUpDown6TissuesDerivedCells_BloodCellsUp,
          flat_allGrSignatures_18.11.2025$minusGlobalUpDown6TissuesDerivedCells_BloodCellsDown
)

# neural
flat_allGrSignatures_18.11.2025$minusGlobalUpDown6TissuesDerivedCells_NeuralCellsUp %>% length()
flat_allGrSignatures_18.11.2025$minusGlobalUpDown6TissuesDerivedCells_NeuralCellsDown %>% length()

intersect(flat_allGrSignatures_18.11.2025$minusGlobalUpDown6TissuesDerivedCells_NeuralCellsUp,
          flat_allGrSignatures_18.11.2025$minusGlobalUpDown6TissuesDerivedCells_NeuralCellsDown
)

# lung
flat_allGrSignatures_18.11.2025$minusGlobalUpDown6TissuesDerivedCells_LungCellsUp %>% length()
flat_allGrSignatures_18.11.2025$minusGlobalUpDown6TissuesDerivedCells_LungCellsDown %>% length()

intersect(flat_allGrSignatures_18.11.2025$minusGlobalUpDown6TissuesDerivedCells_LungCellsUp,
          flat_allGrSignatures_18.11.2025$minusGlobalUpDown6TissuesDerivedCells_LungCellsDown
)



# ##############################################################################
# ---- save to xlsx ----
# ##############################################################################

library(openxlsx)
library(dplyr)

# ============================================================
#  Formatting style + function to save XLSX with formatting
# ============================================================

header_style <- createStyle(textDecoration = "bold", halign = "center")

save_formatted_xlsx <- function(df, path) {
  wb <- createWorkbook()
  addWorksheet(wb, "genes")
  writeData(wb, sheet = 1, x = df, headerStyle = header_style)
  freezePane(wb, sheet = 1, firstActiveRow = 2)   # lock header row
  setColWidths(wb, sheet = 1, cols = 1:ncol(df), widths = "auto")
  saveWorkbook(wb, path, overwrite = TRUE)
  cat("Zapisano:", path, "\n")
}

# ======================================================================
# =============================  MINUS 4  ===============================
# ======================================================================

df_minus4_list <- list(
  globalUp        = flat_allGrSignatures_18.11.2025$global_GR_genes_globalUp4TissuesDerivedCells,
  globalDown      = flat_allGrSignatures_18.11.2025$global_GR_genes_globalDown4TissuesDerivedCells,
  NeuralCellsUp   = flat_allGrSignatures_18.11.2025$minusGlobalUpDown4TissuesDerivedCells_NeuralCellsUp,
  NeuralCellsDown = flat_allGrSignatures_18.11.2025$minusGlobalUpDown4TissuesDerivedCells_NeuralCellsDown,
  BloodCellsUp    = flat_allGrSignatures_18.11.2025$minusGlobalUpDown4TissuesDerivedCells_BloodCellsUp,
  BloodCellsDown  = flat_allGrSignatures_18.11.2025$minusGlobalUpDown4TissuesDerivedCells_BloodCellsDown,
  LungCellsUp     = flat_allGrSignatures_18.11.2025$minusGlobalUpDown4TissuesDerivedCells_LungCellsUp,
  LungCellsDown   = flat_allGrSignatures_18.11.2025$minusGlobalUpDown4TissuesDerivedCells_LungCellsDown
)

max_len4 <- max(sapply(df_minus4_list, length))
df_minus4 <- lapply(df_minus4_list, function(x) c(x, rep(NA, max_len4 - length(x)))) %>% as.data.frame()

path4 <- "/home/mateusz/projects/ifpan-michkor-gr/results_v2/positive_validation_grSignatures_18.11.2025/flat_geneLists_minus4_global_tissue.xlsx"
save_formatted_xlsx(df_minus4, path4)

# ======================================================================
# =============================  MINUS 5  ===============================
# ======================================================================

df_minus5_list <- list(
  globalUp        = flat_allGrSignatures_18.11.2025$global_GR_genes_globalUp5TissuesDerivedCells,
  globalDown      = flat_allGrSignatures_18.11.2025$global_GR_genes_globalDown5TissuesDerivedCells,
  NeuralCellsUp   = flat_allGrSignatures_18.11.2025$minusGlobalUpDown5TissuesDerivedCells_NeuralCellsUp,
  NeuralCellsDown = flat_allGrSignatures_18.11.2025$minusGlobalUpDown5TissuesDerivedCells_NeuralCellsDown,
  BloodCellsUp    = flat_allGrSignatures_18.11.2025$minusGlobalUpDown5TissuesDerivedCells_BloodCellsUp,
  BloodCellsDown  = flat_allGrSignatures_18.11.2025$minusGlobalUpDown5TissuesDerivedCells_BloodCellsDown,
  LungCellsUp     = flat_allGrSignatures_18.11.2025$minusGlobalUpDown5TissuesDerivedCells_LungCellsUp,
  LungCellsDown   = flat_allGrSignatures_18.11.2025$minusGlobalUpDown5TissuesDerivedCells_LungCellsDown
)

max_len5 <- max(sapply(df_minus5_list, length))
df_minus5 <- lapply(df_minus5_list, function(x) c(x, rep(NA, max_len5 - length(x)))) %>% as.data.frame()

path5 <- "/home/mateusz/projects/ifpan-michkor-gr/results_v2/positive_validation_grSignatures_18.11.2025/flat_geneLists_minus5_global_tissue.xlsx"
save_formatted_xlsx(df_minus5, path5)

# ======================================================================
# =============================  MINUS 6  ===============================
# ======================================================================

df_minus6_list <- list(
  globalUp        = flat_allGrSignatures_18.11.2025$global_GR_genes_globalUp6TissuesDerivedCells,
  globalDown      = flat_allGrSignatures_18.11.2025$global_GR_genes_globalDown6TissuesDerivedCells,
  NeuralCellsUp   = flat_allGrSignatures_18.11.2025$minusGlobalUpDown6TissuesDerivedCells_NeuralCellsUp,
  NeuralCellsDown = flat_allGrSignatures_18.11.2025$minusGlobalUpDown6TissuesDerivedCells_NeuralCellsDown,
  BloodCellsUp    = flat_allGrSignatures_18.11.2025$minusGlobalUpDown6TissuesDerivedCells_BloodCellsUp,
  BloodCellsDown  = flat_allGrSignatures_18.11.2025$minusGlobalUpDown6TissuesDerivedCells_BloodCellsDown,
  LungCellsUp     = flat_allGrSignatures_18.11.2025$minusGlobalUpDown6TissuesDerivedCells_LungCellsUp,
  LungCellsDown   = flat_allGrSignatures_18.11.2025$minusGlobalUpDown6TissuesDerivedCells_LungCellsDown
)

max_len6 <- max(sapply(df_minus6_list, length))
df_minus6 <- lapply(df_minus6_list, function(x) c(x, rep(NA, max_len6 - length(x)))) %>% as.data.frame()

path6 <- "/home/mateusz/projects/ifpan-michkor-gr/results_v2/positive_validation_grSignatures_18.11.2025/flat_geneLists_minus6_global_tissue.xlsx"
save_formatted_xlsx(df_minus6, path6)

# ======================================================================
# =============== CLEANUP: delete helper variables from environment =====
# ======================================================================

rm(list = ls(pattern = "df_minus|max_len|path[456]|_list|header_style|save_formatted_xlsx"))

cat("\n🧹 Usunięto zmienne pomocnicze z environment.\n")

