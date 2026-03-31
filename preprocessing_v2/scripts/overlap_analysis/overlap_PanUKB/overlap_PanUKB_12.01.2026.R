list.files(
  "data/prs-models-pan-biobank-uk/",
  pattern = "1e-08\\.yml$",
  full.names = TRUE
) %>%
  # Set the names of the list elements to the file names without the extension
  set_names(map(., ~ basename(.x) %>% tools::file_path_sans_ext())) %>%
  # Read YAML content from each file
  map(read_yaml) %>%
  # Extract genes information and name the list elements after the files
  map(~ .x$description$genes) %>% lapply(., unique) -> phenotypes_PanUkBiobank


phenotypes_PanUkBiobank %>% lapply(., length)


keep(
  phenotypes_PanUkBiobank,
  ~ length(.x) >= 10
) -> gene_list


PanUKB_GrSignatures_overlapChi2 <- run_full_overlap_analysis(
  gene_lists = c(flat_allGrSignatures_31.10.2025[sig_names], 
                 gene_list
  ),
  total_genes = hgnc_symbols_vector_v110,
  cols_to_filter = sig_names,
  rows_to_filter = names(gene_list),
  plot_title_or = "",
  triangle_mode = "full",
  fdr_threshold = 1,
  data_type = "original_data",
  verbose = FALSE,
  palette_or = c("#c6d3e3", "white", "darkred"),
  text_contrast_range_or = c(-30, 4.9),
)


PanUKB_GrSignatures_overlapChi2$processed$original_data$df %>% 
  select(-c(sig_1, ._key)) %>% 
  filter(fdr_value < 0.01) %>% 
  filter(gene_overlap_count  >= 3) %>% 
  filter(log2_odds_ratio > 0) %>% 
  .$Var2 %>% as.character() %>% unique() -> PanUKB_grSignatures_vector_fdr0.01


PanUKB_GrSignatures_overlapChi2$processed$original_data$df %>% 
  select(-c(sig_1, ._key)) %>% 
  filter(fdr_value < 0.01) %>% 
  filter(gene_overlap_count  >= 3) %>% 
  filter(log2_odds_ratio > 0) %>% 
  .$Var1 %>% as.character() %>% unique() -> PanUKB_phenotypes_vector_fdr0.01


heatmap_overlap_log2OR_complex(
  data_list = PanUKB_GrSignatures_overlapChi2$processed,
  data_type = "original_data",
  
  # 🔺 triangle mode pełny
  triangle_mode = "full",
  
  # 🔹 filtrowanie
  rows_to_filter = PanUKB_phenotypes_vector_fdr0.01,
  # cols_to_filter = pgc_grSignatures_vector_p0.05,
  
  p_thresholds = c(0.005, 0.00001),
  
  # 🎨 skala kolorów
  color_scale_range = c(-3, 3),
  text_contrast_range = c(-30, 4.9),
  palette = c("#c6d3e3", "white", "darkred"),
  
  # 🔲 prostokąty istotnościhttp://localhost:8600/graphics/plot_zoom_png?width=1029&height=2556
  color_rects = c("#97C426", "#2F4603"),
  
  # 📊 klastrowanie
  cluster_rows = T,
  cluster_cols = F,
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
  # save_to_svg = "results_v2/overlap/pgc_overlap/figures/heatmap_allGrSignaturesPGC_log2OR.svg",
  svg_width = 10.5, 
  svg_height = 10,
  row_names_width = unit(90, "mm"),
  col_names_height = unit(50, "mm"),
  force_create_directory = TRUE
)



df_list_by_gr <- PanUKB_GrSignatures_overlapChi2$processed$original_data$df %>% 
  select(-c(sig_1, ._key, fdr)) %>% 
  filter(fdr_value < 0.01) %>% 
  filter(gene_overlap_count >= 3) %>% 
  filter(log2_odds_ratio > 0) %>%  
  mutate(Var2 = str_replace_all(Var2, "minusGlobalUpDown5TissuesDerivedCells_", "")) %>% 
  mutate(Var2 = str_replace_all(Var2, "global_GR_genes_globalUp5TissuesDerivedCells", "globalUp")) %>% 
  mutate(Var2 = str_replace_all(Var2, "global_GR_genes_globalDown5TissuesDerivedCells", "globalDown")) %>% 
  rename(grSignature = Var2,
         PanUKBphenotype = Var1) %>% 
  
  group_by(grSignature) %>% 
  group_split()

# nadajemy sensowne nazwy elementom listy
names(df_list_by_gr) <- unique(
  PanUKB_GrSignatures_overlapChi2$processed$original_data$df %>% 
    mutate(Var2 = str_replace_all(Var2, "minusGlobalUpDown5TissuesDerivedCells_", "")) %>% 
    mutate(Var2 = str_replace_all(Var2, "global_GR_genes_globalUp5TissuesDerivedCells", "globalUp")) %>% 
    mutate(Var2 = str_replace_all(Var2, "global_GR_genes_globalDown5TissuesDerivedCells", "globalDown")) %>% 
    pull(Var2)
)



suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(openxlsx)
})

# ============================================================
# 1) Build list of data.frames split by grSignature
#    (based on your pipeline)
# ============================================================
df_list_by_gr <- PanUKB_GrSignatures_overlapChi2$processed$original_data$df %>%
  select(-c(sig_1, ._key, fdr)) %>%
  filter(fdr_value < 0.01) %>%
  filter(gene_overlap_count >= 3) %>%
  filter(log2_odds_ratio > 0) %>%
  mutate(
    Var2 = str_replace_all(Var2, "minusGlobalUpDown5TissuesDerivedCells_", ""),
    Var2 = str_replace_all(Var2, "global_GR_genes_globalUp5TissuesDerivedCells", "globalUp"),
    Var2 = str_replace_all(Var2, "global_GR_genes_globalDown5TissuesDerivedCells", "globalDown")
  ) %>%
  rename(grSignature = Var2,
         PanUKBphenotype = Var1) %>%
  split(.$grSignature)

# ============================================================
# 2) Helpers for Excel output
# ============================================================
sanitize_sheet_name <- function(x) {
  x <- gsub("[:\\\\/\\?\\*\\[\\]]", "_", x) # forbidden in Excel
  x <- trimws(x)
  x <- substr(x, 1, 31)                     # max 31 chars
  if (nchar(x) == 0) x <- "sheet"
  x
}

write_list_to_xlsx_with_formatting <- function(df_list, file, fdr_col = "fdr_value") {
  stopifnot(is.list(df_list), length(df_list) > 0)
  
  wb <- createWorkbook()
  green_row_style <- createStyle(fgFill = "#C6EFCE")  # light green
  
  used_names <- character(0)
  
  for (nm in names(df_list)) {
    df <- df_list[[nm]]
    if (!is.data.frame(df)) next
    
    sheet <- sanitize_sheet_name(nm)
    
    # ensure unique sheet names
    base <- sheet
    k <- 1
    while (sheet %in% used_names) {
      suffix <- paste0("_", k)
      sheet <- substr(paste0(substr(base, 1, 31 - nchar(suffix)), suffix), 1, 31)
      k <- k + 1
    }
    used_names <- c(used_names, sheet)
    
    addWorksheet(wb, sheetName = sheet)
    
    # write
    writeData(wb, sheet = sheet, x = df, withFilter = TRUE)
    
    # freeze first row + first column
    freezePane(wb, sheet = sheet, firstActiveRow = 2, firstActiveCol = 2)
    
    # auto width
    setColWidths(wb, sheet = sheet, cols = 1:ncol(df), widths = "auto")
    
    # green rows where fdr_value < 0.05
    if (fdr_col %in% colnames(df) && nrow(df) > 0) {
      fdr_idx <- match(fdr_col, colnames(df))
      fdr_letter <- int2col(fdr_idx)
      
      start_row <- 2
      end_row   <- nrow(df) + 1
      start_col <- 1
      end_col   <- ncol(df)
      
      rule_formula <- paste0("$", fdr_letter, start_row, "<0.05")
      
      conditionalFormatting(
        wb, sheet = sheet,
        cols = start_col:end_col,
        rows = start_row:end_row,
        rule = rule_formula,
        style = green_row_style,
        type = "expression"
      )
    }
  }
  
  saveWorkbook(wb, file = file, overwrite = TRUE)
  invisible(file)
}

# ============================================================
# 3) Save here:
# ============================================================
out_dir  <- "/home/mateusz/projects/ifpan-michkor-gr/results_v2/overlap/PanUKB_ovelap"
out_file <- file.path(out_dir, "PanUKBBiobank_GrSignatures_12.01.2026.xlsx")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

write_list_to_xlsx_with_formatting(
  df_list = df_list_by_gr,
  file    = out_file,
  fdr_col = "fdr_value"
)

out_file

