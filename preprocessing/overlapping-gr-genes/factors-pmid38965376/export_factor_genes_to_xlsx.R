# ##############################################################################
# ---- uses data ----
# ##############################################################################
factors_rsidGenes_P1e2locus50kb %>% head


# ##############################################################################
# ---- preprocessing data ----
# ##############################################################################

toCamelCase <- function(snake_str) {
  parts <- strsplit(snake_str, "_")[[1]]
  paste0(parts[1], paste0(toupper(substring(parts[-1], 1,1)), substring(parts[-1], 2)), collapse = "")
}

# 1. Prepare ranked list with categories
preprocessing_list_ranked <- split(factors_rsidGenes_P1e2geneCenter50kb, factors_rsidGenes_P1e2geneCenter50kb$factor_name) %>% 
  lapply(function(x) {
    x %>% 
      dplyr::filter(pvalue < 0.0001) %>% 
      dplyr::arrange(pvalue) %>% 
      dplyr::mutate(
        rank = dplyr::row_number(),
        rank_category = dplyr::case_when(
          rank <= 10  ~ "top10",
          rank <= 50  ~ "top50",
          rank <= 100 ~ "top100",
          rank <= 200 ~ "top200",
          TRUE        ~ "all"
        )
      )
  })

# 2. Prepare sheet names: convert to camelCase, trim, and make unique
sheet_names_raw <- names(preprocessing_list_ranked)
sheet_names_camel <- sapply(sheet_names_raw, toCamelCase)
sheet_names <- substr(sheet_names_camel, 1, 31)
sheet_names <- make.unique(sheet_names, sep = "_")

# 3. Create workbook and write sheets
wb <- openxlsx::createWorkbook()

# Proste jasnozielone kolory (ciemniejszy dla top200, jaśniejszy dla top10)
light_green_colors <- c(
  "top10" = "#9ccc65",  # light-green lighten-1
  "top50" = "#aed581",  # light-green lighten-2
  "top100"  = "#c5e1a5",  # light-green lighten-3
  "top200"  = "#dcedc8"   # light-green lighten-4
)

for (i in seq_along(preprocessing_list_ranked)) {
  openxlsx::addWorksheet(wb, sheet_names[i])
  openxlsx::writeData(wb, sheet = sheet_names[i], x = preprocessing_list_ranked[[i]])
  openxlsx::freezePane(wb, sheet = sheet_names[i], firstRow = TRUE)
  openxlsx::addFilter(wb, sheet = sheet_names[i], row = 1, cols = 1:ncol(preprocessing_list_ranked[[i]]))
  
  rank_cat_col <- which(colnames(preprocessing_list_ranked[[i]]) == "rank_category")
  total_rows <- nrow(preprocessing_list_ranked[[i]]) + 1
  
  first_col <- 1
  last_col <- ncol(preprocessing_list_ranked[[i]])
  
  addColorRule <- function(sheet, value, color) {
    openxlsx::conditionalFormatting(
      wb, sheet = sheet, cols = first_col:last_col, rows = 2:(total_rows),
      rule = paste0('$', openxlsx::int2col(rank_cat_col), '2="', value, '"'),
      style = openxlsx::createStyle(bgFill = color),
      type = "expression"
    )
  }
  
  addColorRule(sheet_names[i], "top10",  light_green_colors["top10"])
  addColorRule(sheet_names[i], "top50",  light_green_colors["top50"])
  addColorRule(sheet_names[i], "top100", light_green_colors["top100"])
  addColorRule(sheet_names[i], "top200", light_green_colors["top200"])
}

# 4. Save workbook with descriptive filename
openxlsx::saveWorkbook(wb, "data/factors-pmid38965376/factorsGeneList_P1e4geneCenter50kb.xlsx", overwrite = TRUE)

# 5. Clean up
rm(preprocessing_list_ranked, wb, sheet_names_raw, sheet_names_camel, sheet_names)
