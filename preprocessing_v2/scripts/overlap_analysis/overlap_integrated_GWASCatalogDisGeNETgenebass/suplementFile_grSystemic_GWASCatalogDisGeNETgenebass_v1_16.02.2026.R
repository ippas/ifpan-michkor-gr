GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$df %>% 
  filter(!(Var1 %in% c("DisGeNET_NA", "genebass_NA", "GWASCatalog_NA",
                       "genebass_F9x_Date_F99_first_reported_(mental_disorder,_not_otherwise_specified)"))) %>% 
  rowwise() %>%
  mutate(
    matched = list(str_extract_all(Var1, paste(patterns, collapse="|"))[[1]]),
    n_unique_patterns = n_distinct(matched)
  ) %>%
  ungroup() %>%
  filter(n_unique_patterns < 2) %>%
  filter(n_unique_patterns != 0) %>% 
  mutate(
    icd10_category = map_chr(matched, ~ {
      if (length(.x) == 1) .x else NA_character_
    })
  ) %>%
  mutate(
    Var2 = recode(Var2,
                  "global_GR_genes_globalUp5TissuesDerivedCells"   = "systemicUp",
                  "global_GR_genes_globalDown5TissuesDerivedCells" = "systemicDown"
    )
  ) %>% 
  select(-c(fdr, fdr_value, matched, n_unique_patterns)) %>% 

  separate(Var1,
           into = c("source", "phenotype"),
           sep = "_",
           extra = "merge") %>% 
  mutate(
    phenotype = str_remove(phenotype, "^([A-Z][0-9]x_)+")
  ) %>% as.data.frame() %>% 
  rename(grSignature = "Var2") %>% 
  rename(grSignature_nGenes = "Var2_n_genes") %>% 
  rename(phenotype_nGenes = "Var1_n_genes") %>% 
  mutate(minusLog10_pvalue = -log10(p_value)) %>% 
  mutate(log10_overlapPlusOne = log10(gene_overlap_count + 1)) %>% 
  mutate(combine_score = minusLog10_pvalue * log10_overlapPlusOne) %>% 
  mutate(
    chi2_residual = (gene_overlap_count - expected_overlap) / sqrt(expected_overlap)
  ) %>% 
  select(
    icd10_category,
    phenotype,
    source,
    
    phenotype_nGenes,
    grSignature,
    grSignature_nGenes,
    expected_overlap,
    
    gene_overlap_count,
    
    odds_ratio,
    
    chi2_residual,
    chi2,
    p_value,
    
    minusLog10_pvalue,
    log10_overlapPlusOne,
    combine_score,
    overlap_genes
  ) %>% 
  rename(observed_overlap = "gene_overlap_count") %>%
  select(-c(chi2_residual, minusLog10_pvalue, log10_overlapPlusOne)) %>% 
  mutate(
    grSignature = factor(grSignature,
                         levels = c("systemicUp", "systemicDown"))
  ) %>%
  arrange(grSignature, desc(combine_score)) %>% 
  split(.$icd10_category) -> grSystemic_GWASCatalogDisGeNETgenebass_list
  

# ##############################################################################
# ---- save to file ----
# ##############################################################################
# -------------------------
# input
# -------------------------
df_list <- grSystemic_GWASCatalogDisGeNETgenebass_list

out_dir  <- "/home/mateusz/projects/ifpan-michkor-gr/results_v2/overlap/overlap_GWASCatalogDisGeNETgenebass/tables"
out_file <- file.path(out_dir, "overlap_grSystemic_GWASCatalogDisGeNETgenebassICD10F_v1_16.02.2026.xlsx")

# kolumny do formatowania (3 miejsca po przecinku)
fmt_3dec_cols <- c(
  "minusLog10_pvalue",
  "expected_overlap",
  "odds_ratio",
  "chi2_residual",
  "chi2",
  "combine_score",
  "log10_overlapPlusOne"
)

# -------------------------
# workbook
# -------------------------
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
wb <- createWorkbook()

style_3dec <- createStyle(numFmt = "0.000")
style_sci3 <- createStyle(numFmt = "0.000E+00")
style_green <- createStyle(fgFill = "#C6EFCE")  # jasnozielony (Excel-like)

# helper: bezpieczne nazwy arkuszy (<= 31 znaków)
safe_sheet_name <- function(x) {
  x <- gsub("[:\\\\/?*\\[\\]]", "_", x)  # znaki niedozwolone
  x <- substr(x, 1, 31)
  x
}

# -------------------------
# write each df to separate sheet
# -------------------------
sheet_names <- names(df_list)
if (is.null(sheet_names) || any(sheet_names == "")) {
  sheet_names <- paste0("sheet_", seq_along(df_list))
}

sheet_names <- vapply(sheet_names, safe_sheet_name, character(1))

for (i in seq_along(df_list)) {
  df <- df_list[[i]]
  sh <- sheet_names[i]
  
  addWorksheet(wb, sh)
  
  # zapis danych
  writeData(wb, sh, df, withFilter = TRUE)
  
  # zamrożenie 1. wiersza
  freezePane(wb, sh, firstRow = TRUE)
  
  # auto szerokości (opcjonalnie; możesz usunąć jak pliki będą zbyt ciężkie)
  setColWidths(wb, sh, cols = 1:ncol(df), widths = "auto")
  
  # --- formaty liczb ---
  # 3 miejsca po przecinku
  cols_3dec <- intersect(fmt_3dec_cols, names(df))
  if (length(cols_3dec) > 0) {
    addStyle(
      wb, sh, style_3dec,
      rows = 2:(nrow(df) + 1),
      cols = match(cols_3dec, names(df)),
      gridExpand = TRUE,
      stack = TRUE
    )
  }
  
  # p_value: zapis naukowy, 3 miejsca po przecinku
  if ("p_value" %in% names(df)) {
    pcol <- match("p_value", names(df))
    addStyle(
      wb, sh, style_sci3,
      rows = 2:(nrow(df) + 1),
      cols = pcol,
      gridExpand = TRUE,
      stack = TRUE
    )
  }
  
  # --- podświetlenie: p_value < 0.05 ORAZ observed_overlap >= 3 ---
  # podświetlamy komórki p_value (jeśli wolisz cały wiersz, daj znać)
  if (all(c("p_value", "observed_overlap") %in% names(df)) && nrow(df) > 0) {
    pcol <- match("p_value", names(df))
    ocol <- match("observed_overlap", names(df))
    
    # Excel formula; 2 = pierwszy wiersz danych (po nagłówku)
    p_letter <- int2col(pcol)
    o_letter <- int2col(ocol)
    rule <- paste0("AND($", p_letter, "2<0.05,$", o_letter, "2>=3)")
    
    conditionalFormatting(
      wb, sh,
      cols = 1:ncol(df),          # <- cały wiersz
      rows = 2:(nrow(df) + 1),
      rule = rule,
      style = style_green,
      type = "expression"
    )
  }
}

saveWorkbook(wb, out_file, overwrite = TRUE)

out_file


# ##############################################################################
grSystemic_GWASCatalogDisGeNETgenebass_list$F8x %>% 
  filter(p_value < 0.05, observed_overlap >= 3, odds_ratio > 1) %>% .$combine_score %>%   paste(collapse = "|")

grSystemic_GWASCatalogDisGeNETgenebass_list$F8x %>% 
  filter(p_value < 0.05, observed_overlap >= 3, odds_ratio > 1) %>% .$combine_score %>% mean

grSystemic_GWASCatalogDisGeNETgenebass_list$F8x %>% 
  filter(p_value < 0.05, observed_overlap >= 3, odds_ratio > 1) %>% pull(overlap_genes) %>% 
  strsplit(",") %>%
  unlist %>% 
  unique %>% 
  sort %>%
  paste(collapse = "|")

grSystemic_GWASCatalogDisGeNETgenebass_list %>% 
  bind_rows() %>% 
  filter(p_value < 0.05) %>% 
  # filter(icd10_category %in% c("F3x")) %>% 
  select(grSignature, combine_score) %>% 
  t.test(combine_score ~ grSignature, data = .)
  
grSystemic_GWASCatalogDisGeNETgenebass_list %>% 
  bind_rows() %>% 
  filter(p_value < 0.05 & observed_overlap > 2) %>% 
  select(grSignature, combine_score) %>% 
  ggplot(aes(x = grSignature, y = combine_score, fill = grSignature)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  theme_classic() +
  labs(
    x = "GR signature",
    y = "Combined score",
    title = "Combined score by GR signature (p < 0.05)"
  ) +
  theme(legend.position = "none")


# ##############################################################################
grSystemic_GWASCatalogDisGeNETgenebass_list$F3x %>% 
  filter(p_value < 0.05) %>% 
  filter(observed_overlap >= 3)
