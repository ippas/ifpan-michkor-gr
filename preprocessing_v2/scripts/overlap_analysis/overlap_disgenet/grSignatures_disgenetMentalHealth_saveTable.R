# ======================================
# PREPARE DATA
# ======================================
df_full <- disgenetMentalHealth_GrSignatures_overlapChi2$processed$original_data$df %>% 
  select(-c(sig_1, fdr_value, fdr, ._key)) %>% 
  mutate(Var2 = str_replace_all(Var2, "minusGlobalUpDown5TissuesDerivedCells_", "")) %>% 
  mutate(Var2 = ifelse(Var2 == "global_GR_genes_globalUp5TissuesDerivedCells", "globalUp", Var2)) %>% 
  mutate(Var2 = ifelse(Var2 == "global_GR_genes_globalDown5TissuesDerivedCells", "globalDown", Var2)) %>% 
  rename(Var2 = "grSignature") %>% 
  rename(Var1 = "GWASfile")

df_significant <- df_full %>%
  filter(
    gene_overlap_count > 2,
    p_value < 0.05,
    log2_odds_ratio > 0
  )

# ======================================
# CREATE EXCEL WORKBOOK
# ======================================
wb <- createWorkbook()

addWorksheet(wb, "significant_results")
writeData(wb, "significant_results", df_significant)
freezePane(wb, "significant_results", firstActiveRow = 2, firstActiveCol = 1)

addWorksheet(wb, "full_results")
writeData(wb, "full_results", df_full)
freezePane(wb, "full_results", firstActiveRow = 2, firstActiveCol = 1)

# Styles
style_numeric <- createStyle(numFmt = "0.000")
style_sci <- createStyle(numFmt = "0.000E+00")

numeric_cols <- c("odds_ratio", "log2_odds_ratio", "chi2")
pvalue_col <- "p_value"

# Apply formatting & widths
for(sheet in c("significant_results", "full_results")) {
  for(col in numeric_cols){
    if(col %in% names(df_full)){
      addStyle(wb, sheet, style_numeric,
               cols = which(names(df_full) == col),
               rows = 2:(nrow(df_full)+1),
               gridExpand = TRUE)
    }
  }
  
  addStyle(wb, sheet, style_sci,
           cols = which(names(df_full) == pvalue_col),
           rows = 2:(nrow(df_full)+1),
           gridExpand = TRUE)
  
  setColWidths(wb, sheet, cols = 1:ncol(df_full), widths = "auto")
}

# SAVE FILE
output_path <- "/home/mateusz/projects/ifpan-michkor-gr/results_v2/overlap/disgenet_overlap/tables/overlap_grSignaturesDisGeNET_chi2_p0.05minOverlapGenes3-16.11.2025.xlsx"

saveWorkbook(wb, output_path, overwrite = TRUE)

message("💾 Saved Excel file to: ", output_path)

