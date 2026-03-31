# ##############################################################################
# ---- data ----
# ##############################################################################
sig_names <- c(
  "minusGlobalUpDown5TissuesDerivedCells_BloodCellsUp",
  "minusGlobalUpDown5TissuesDerivedCells_LungCellsUp",
  "minusGlobalUpDown5TissuesDerivedCells_NeuralCellsUp",
  "minusGlobalUpDown5TissuesDerivedCells_BloodCellsDown",
  "minusGlobalUpDown5TissuesDerivedCells_LungCellsDown",
  "minusGlobalUpDown5TissuesDerivedCells_NeuralCellsDown",
  "global_GR_genes_globalUp5TissuesDerivedCells",
  "global_GR_genes_globalDown5TissuesDerivedCells"
)

flat_allGrSignatures_31.10.2025[sig_names]


# ##############################################################################
# ---- run enrichr ----
# ##############################################################################

enrichr_GWASCatalog_grSystemicTissues <- run_enrichr_multi(gene_lists = flat_allGrSignatures_31.10.2025[sig_names],
                  database = "GWAS_Catalog_2025",
                  min_overlap_genes = 0,
                  fdr_threshold = 1)


new_names <- c(
  bloodUp     = "minusGlobalUpDown5TissuesDerivedCells_BloodCellsUp",
  bloodDown   = "minusGlobalUpDown5TissuesDerivedCells_BloodCellsDown",
  neuralUp    = "minusGlobalUpDown5TissuesDerivedCells_NeuralCellsUp",
  neuralDown  = "minusGlobalUpDown5TissuesDerivedCells_NeuralCellsDown",
  lungUp      = "minusGlobalUpDown5TissuesDerivedCells_LungCellsUp",
  lungDown    = "minusGlobalUpDown5TissuesDerivedCells_LungCellsDown",
  systemicUp  = "global_GR_genes_globalUp5TissuesDerivedCells",
  systemicDown= "global_GR_genes_globalDown5TissuesDerivedCells"
)

names(enrichr_GWASCatalog_grSystemicTissues$enrichr$raw) <-
  names(new_names)[match(
    names(enrichr_GWASCatalog_grSystemicTissues$enrichr$raw),
    new_names
  )]

enrichr_GWASCatalog_grSystemicTissues$enrichr$raw %>% 
  lapply(., function(x){
    x %>% 
      filter(n_genes >= 3) %>% 
      filter(pvalue < 0.05)
  })



# -----------------------------------------------------------------------------
# Desired sheet order
# -----------------------------------------------------------------------------
library(dplyr)
library(openxlsx)

# -----------------------------------------------------------------------------
# Desired sheet order
# -----------------------------------------------------------------------------
sheet_order <- c(
  "systemicUp",
  "systemicDown",
  "neuralUp",
  "neuralDown",
  "bloodUp",
  "bloodDown",
  "lungUp",
  "lungDown"
)

# -----------------------------------------------------------------------------
# Filter enrichr results
# -----------------------------------------------------------------------------
filtered_results <- enrichr_GWASCatalog_grSystemicTissues$enrichr$raw %>%
  lapply(function(x){
    x %>%
      filter(n_genes >= 3) %>%
      filter(pvalue < 0.05) %>%
      select(Term, n_genes, everything())   # move n_genes to column 2
  })

filtered_results <- filtered_results[sheet_order]

# -----------------------------------------------------------------------------
# Create workbook
# -----------------------------------------------------------------------------
wb <- createWorkbook()

# style for significant rows
sig_style <- createStyle(
  fgFill = "#e9edc9",
  border = "TopBottomLeftRight",
  borderColour = "#D9D9D9"
)

for(sheet in sheet_order){
  
  df <- filtered_results[[sheet]]
  
  addWorksheet(
    wb,
    sheetName = sheet,
    gridLines = TRUE
  )
  
  writeData(
    wb,
    sheet,
    df,
    withFilter = TRUE
  )
  
  # freeze header row
  freezePane(
    wb,
    sheet,
    firstRow = TRUE
  )
  
  # ----------------------------------------------------------------------------
  # numeric formatting
  # ----------------------------------------------------------------------------
  
  if("Odds.Ratio" %in% colnames(df)){
    addStyle(
      wb,
      sheet,
      createStyle(numFmt = "0.000"),
      rows = 2:(nrow(df)+1),
      cols = which(colnames(df) == "Odds.Ratio"),
      gridExpand = TRUE
    )
  }
  
  if("Combined.Score" %in% colnames(df)){
    addStyle(
      wb,
      sheet,
      createStyle(numFmt = "0.000"),
      rows = 2:(nrow(df)+1),
      cols = which(colnames(df) == "Combined.Score"),
      gridExpand = TRUE
    )
  }
  
  if("pvalue" %in% colnames(df)){
    addStyle(
      wb,
      sheet,
      createStyle(numFmt = "0.000E+00"),
      rows = 2:(nrow(df)+1),
      cols = which(colnames(df) == "pvalue"),
      gridExpand = TRUE
    )
  }
  
  if("FDR" %in% colnames(df)){
    addStyle(
      wb,
      sheet,
      createStyle(numFmt = "0.000E+00"),
      rows = 2:(nrow(df)+1),
      cols = which(colnames(df) == "FDR"),
      gridExpand = TRUE
    )
    
    # highlight rows where FDR < 0.05
    conditionalFormatting(
      wb,
      sheet,
      cols = 1:ncol(df),
      rows = 2:(nrow(df)+1),
      rule = paste0("$", LETTERS[which(colnames(df)=="FDR")], "2<0.05"),
      style = sig_style
    )
  }
  
  setColWidths(
    wb,
    sheet,
    cols = 1:ncol(df),
    widths = "auto"
  )
}

# -----------------------------------------------------------------------------
# Save workbook
# -----------------------------------------------------------------------------
output_path <- "/home/mateusz/projects/ifpan-michkor-gr/results_v2/enrichr/enrichr_GWASCatalog_grSystemicTissues_13.03.2026/enrichr_GWASCatalog_grSystemicTissues_p0.05genes3__v1_13.03.2023.xlsx"

saveWorkbook(
  wb,
  output_path,
  overwrite = TRUE
)

