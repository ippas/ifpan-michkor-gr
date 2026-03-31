# ============================================================
# Export overlap table to Excel (supplement) — FULL SCRIPT
# - mapped_icd10_range relabeled:
#     F0x -> F00-F09 ... F8x -> F80-F89, F9x -> F90-F98
# - split into sheets by mapped_icd10_range (ordered)
# - within each sheet:
#     * order by grSignature:
#       systemicUp, systemicDown, neuralUp, neuralDown, bloodUp, bloodDown, lungUp, lungDown
#     * then sort by combine_score (DESC)
# - freeze first row + filter
# - highlight rows where p_value < 0.05 & observed_overlap >= 3 with #C9CBA3
# - add thin borders so grid is visible even on colored rows
# - number formatting:
#     * expected_overlap, combine_score -> 3 decimals
#     * chi2, p_value, odds_ratio -> scientific, 3 decimals
# ============================================================

library(dplyr)
library(stringr)
library(openxlsx)

# ---- INPUT ----
df_in <- GWASCatalogDisGeNETgenebass_GrTissueSystemic_overlapChi2$processed$original_data$df

# ---- helper: relabel ICD10 range ----
icd10_relabel <- function(x) {
  d <- str_match(x, "^F([0-9])x$")[, 2]
  ifelse(
    is.na(d),
    NA_character_,
    ifelse(d == "9", "F90-F98", sprintf("F%s0-F%s9", d, d))
  )
}

# ---- desired order of sheets ----
icd10_sheet_order <- c(
  "F00-F09","F10-F19","F20-F29","F30-F39","F40-F49",
  "F50-F59","F60-F69","F70-F79","F80-F89","F90-F98"
)

# ---- desired order of signatures within each sheet ----
grSignature_order <- c(
  "systemicUp","systemicDown",
  "neuralUp","neuralDown",
  "bloodUp","bloodDown",
  "lungUp","lungDown"
)

# ---- BUILD TABLE ----
df_final <- df_in %>%
  mutate(
    log10_pvalue      = -log10(p_value),
    log10_geneOverlap = log10(gene_overlap_count + 1),
    combine_score     = log10_pvalue * log10_geneOverlap
  ) %>%
  filter(!(Var1 %in% c(
    "DisGeNET_NA", "genebass_NA", "GWASCatalog_NA",
    "genebass_F9x_Date_F99_first_reported_(mental_disorder,_not_otherwise_specified)"
  ))) %>%
  # extract coarse ICD10 bucket and relabel to ranges
  mutate(
    icd10_bucket       = str_extract(Var1, "F[0-9]x"),
    mapped_icd10_range = icd10_relabel(icd10_bucket)
  ) %>%
  filter(!is.na(mapped_icd10_range)) %>%
  select(-icd10_bucket) %>%
  # drop unused cols safely (if present)
  select(-any_of(c("log2_odds_ratio", "fdr_value", "fdr"))) %>%
  rename(
    grSignature        = Var2,
    phenotype          = Var1,
    observed_overlap   = gene_overlap_count,
    grSignature_nGenes = Var2_n_genes,
    phenotype_nGenes   = Var1_n_genes
  ) %>%
  mutate(
    grSignature = case_when(
      str_detect(grSignature, "BloodCellsUp")                  ~ "bloodUp",
      str_detect(grSignature, "BloodCellsDown")                ~ "bloodDown",
      str_detect(grSignature, "LungCellsUp")                   ~ "lungUp",
      str_detect(grSignature, "LungCellsDown")                 ~ "lungDown",
      str_detect(grSignature, "NeuralCellsUp")                 ~ "neuralUp",
      str_detect(grSignature, "NeuralCellsDown")               ~ "neuralDown",
      str_detect(grSignature, "global_GR_genes_globalUp")      ~ "systemicUp",
      str_detect(grSignature, "global_GR_genes_globalDown")    ~ "systemicDown",
      TRUE ~ as.character(grSignature)
    ),
    regulation = case_when(
      str_detect(grSignature, "Up")   ~ "up",
      str_detect(grSignature, "Down") ~ "down",
      TRUE ~ NA_character_
    ),
    grSignature        = factor(grSignature, levels = grSignature_order),
    mapped_icd10_range = factor(mapped_icd10_range, levels = icd10_sheet_order),
    source             = str_extract(phenotype, "^[^_]+"),
    phenotype          = str_remove(phenotype, "^[^_]+_F[0-9]x_")
  ) %>%
  select(
    mapped_icd10_range,
    phenotype,
    source,
    phenotype_nGenes,
    grSignature,
    regulation,
    grSignature_nGenes,
    expected_overlap,
    observed_overlap,
    odds_ratio,
    chi2,
    p_value,
    combine_score,
    overlap_genes
  )

# ---- OUTPUT ----
out_xlsx <- "/home/mateusz/projects/ifpan-michkor-gr/results_v2/overlap/overlap_grTissues_GWASCatalogDisGeNETgenebass/tables/supplementFile_overlap_grSystemicTissues_04.03.2026.xlsx"

# ---- STYLES ----
wb <- createWorkbook()

# highlight significant rows
highlight_style <- createStyle(fgFill = "#C9CBA3")

# thin borders so lines are visible even on filled cells
grid_border_style <- createStyle(
  border = c("top", "bottom", "left", "right"),
  borderStyle = "thin",
  borderColour = "#D9D9D9"
)

# number formats
fmt_3dec <- createStyle(numFmt = "0.000")
fmt_sci3 <- createStyle(numFmt = "0.000E+00")

# ---- WRITE SHEETS ----
for (lvl in icd10_sheet_order) {
  
  df_sheet <- df_final %>%
    filter(mapped_icd10_range == lvl) %>%
    arrange(grSignature, desc(combine_score)) %>%   # order by signature then combine_score
    mutate(grSignature = as.character(grSignature)) # nicer display in Excel
  
  if (nrow(df_sheet) == 0) next
  
  addWorksheet(wb, sheetName = lvl)
  
  writeData(wb, sheet = lvl, x = df_sheet, withFilter = TRUE)
  
  freezePane(wb, sheet = lvl, firstRow = TRUE)
  setColWidths(wb, sheet = lvl, cols = 1:ncol(df_sheet), widths = "auto")
  
  # ranges (include header row)
  all_rows  <- 1:(nrow(df_sheet) + 1)
  all_cols  <- 1:ncol(df_sheet)
  data_rows <- 2:(nrow(df_sheet) + 1)
  
  # 1) borders on whole table (header + data)
  addStyle(
    wb, sheet = lvl, style = grid_border_style,
    rows = all_rows, cols = all_cols,
    gridExpand = TRUE, stack = TRUE
  )
  
  # 2) highlight significant rows (keeps borders)
  highlight_rows <- which(df_sheet$p_value < 0.05 & df_sheet$observed_overlap >= 3)
  if (length(highlight_rows) > 0) {
    addStyle(
      wb, sheet = lvl, style = highlight_style,
      rows = highlight_rows + 1,    # +1 header
      cols = all_cols,
      gridExpand = TRUE, stack = TRUE
    )
  }
  
  # 3) number formatting
  col_expected <- match("expected_overlap", names(df_sheet))
  col_cs       <- match("combine_score", names(df_sheet))
  col_chi2     <- match("chi2", names(df_sheet))
  col_p        <- match("p_value", names(df_sheet))
  col_or       <- match("odds_ratio", names(df_sheet))
  
  addStyle(
    wb, sheet = lvl, style = fmt_3dec,
    rows = data_rows, cols = c(col_expected, col_cs),
    gridExpand = TRUE, stack = TRUE
  )
  
  addStyle(
    wb, sheet = lvl, style = fmt_sci3,
    rows = data_rows, cols = c(col_chi2, col_p, col_or),
    gridExpand = TRUE, stack = TRUE
  )
}

# ---- SAVE ----
saveWorkbook(wb, out_xlsx, overwrite = TRUE)
message("Saved: ", out_xlsx)