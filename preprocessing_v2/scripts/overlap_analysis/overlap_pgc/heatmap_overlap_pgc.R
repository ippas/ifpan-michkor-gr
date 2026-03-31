heatmap_overlap_log2OR_complex(
  data_list = pgcGrSignatures_overlapChi2$processed,
  data_type = "original_data",
  
  # 🔺 triangle mode pełny
  triangle_mode = "full",
  
  # 🔹 filtrowanie
  rows_to_filter = pgc_phenotypes_vector_p0.05,
  # cols_to_filter = pgc_grSignatures_vector_p0.05,
  
  # 🎨 skala kolorów
  color_scale_range = c(-3, 3),
  text_contrast_range = c(-30, 4.9),
  palette = c("#c6d3e3", "white", "darkred"),
  
  # 🔲 prostokąty istotności
  color_rects = c("#97C426", "#2F4603"),
  
  # 📊 klastrowanie
  cluster_rows = T,
  cluster_cols = T,
  show_dendrograms = TRUE,
  rect_lwd = 2.5,
  col_mapper = c(
    "minusGlobalUpDown5TissuesDerivedCells_BloodCellsUp" = "BloodCellsUp",
    "minusGlobalUpDown5TissuesDerivedCells_LungCellsUp" = "LungCellsUp",
    "minusGlobalUpDown5TissuesDerivedCells_NeuralCellsUp" = "NeuralCellsUp",
    "minusGlobalUpDown5TissuesDerivedCells_BloodCellsDown" = "BloodCellsDown",
    "minusGlobalUpDown5TissuesDerivedCells_LungCellsDown" = "LungCellsDown",
    "minusGlobalUpDown5TissuesDerivedCells_NeuralCellsDown" = "NeuralCellsDown",
    "global_GR_genes_globalDown5TissuesDerivedCells" =  "globalDown",
    "global_GR_genes_globalUp5TissuesDerivedCells" =  "globalUp"
  ),
  row_dend_height = unit(20, "mm"),
  col_dend_height = unit(20, "mm"),
  tile_gap = 1,
  # save_to_svg = "results_v2/overlap/pgc_overlap/figures/heatmap_AllGrSignaturesPGC_log2OR.svg",
  svg_width = 10.5, 
  svg_height = 9,
  row_names_width = unit(90, "mm"),
  col_names_height = unit(50, "mm"),
  force_create_directory = TRUE
)



pgcGrSignatures_overlapChi2$processed$original_data$df %>% 
  filter(p_value < 0.05) %>% 
  filter(gene_overlap_count > 2) %>% 
  filter(log2_odds_ratio > 0) %>% 
  select(-._key) %>% 
  filter(Var2 == "global_GR_genes_globalDown5TissuesDerivedCells") %>% 
  select(c(Var1,  gene_overlap_count, overlap_genes))



# ======================================
# PREPARE DATA
# ======================================
df_full <- pgcGrSignatures_overlapChi2$processed$original_data$df %>% 
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

# ---- SIGNIFICANT RESULTS ----
addWorksheet(wb, "significant_results")
writeData(wb, "significant_results", df_significant)
freezePane(wb, "significant_results", firstActiveRow = 2, firstActiveCol = 1)

# ---- FULL RESULTS ----
addWorksheet(wb, "full_results")
writeData(wb, "full_results", df_full)
freezePane(wb, "full_results", firstActiveRow = 2, firstActiveCol = 1)

# ---- STYLES ----
style_numeric <- createStyle(numFmt = "0.000")
style_sci <- createStyle(numFmt = "0.000E+00")

numeric_cols <- c("odds_ratio", "log2_odds_ratio", "chi2")
pvalue_col <- "p_value"

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

# ======================================
# SAVE FILE
# ======================================
output_path <- "/home/mateusz/projects/ifpan-michkor-gr/results_v2/overlap/pgc_overlap/tables/overlap_grSingautresPGC_chi2_p0.05minOverlapGenes3-16.11.2025.xlsx"

saveWorkbook(wb, output_path, overwrite = TRUE)

message("💾 Saved Excel file to: ", output_path)
