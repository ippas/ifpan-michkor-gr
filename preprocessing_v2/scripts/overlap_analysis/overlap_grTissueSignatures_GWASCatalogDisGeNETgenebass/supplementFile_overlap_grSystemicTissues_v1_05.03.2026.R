# ============================================================
# Export overlap table to Excel:
# - compute combine_score
# - map signatures -> (blood/lung/neural/systemic) + up/down
# - extract ICD10 category (F0x..F9x) and then relabel to:
#     F0x -> F00-F09
#     F1x -> F10-F19
#     ...
#     F8x -> F80-F89
#     F9x -> F90-F98   (special case)
# - split into sheets by relabeled ICD10 range (in correct order)
# - freeze first row, add filter, auto column widths
# - save to specified path/name
# ============================================================

library(dplyr)
library(stringr)
library(openxlsx)

# ---- INPUT ----
df_in <-GWASCatalogDisGeNETgenebass_GrTissueSystemic_overlapChi2$processed$original_data$df %>% 
  mutate(
    log10_pvalue = -log10(p_value),
    log10_geneOverlap = log10(gene_overlap_count + 1),
    combine_score = log10_pvalue * log10_geneOverlap
  ) %>%
  filter(!(Var1 %in% c(
    "DisGeNET_NA", "genebass_NA", "GWASCatalog_NA",
    "genebass_F9x_Date_F99_first_reported_(mental_disorder,_not_otherwise_specified)"
  ))) %>%
  rowwise() %>%
  mutate(
    matched = list(str_extract_all(Var1, paste(patterns, collapse="|"))[[1]]),
    n_unique_patterns = n_distinct(matched)
  ) %>%
  ungroup() %>%
  filter(n_unique_patterns == 1) %>% 
  mutate(
    icd10_category = map_chr(matched, ~ .x[[1]])
  ) %>%
  # filter(!is.na(regulation), !is.na(icd10_category)) %>%
  filter(!is.na(icd10_category)) %>% 
  select(-c(log2_odds_ratio, fdr_value, fdr)) %>% 
  rename(
    grSignature        = Var2,
    phenotype          = Var1,
    observed_overlap   = gene_overlap_count,
    grSignature_nGenes = Var2_n_genes,
    phenotype_nGenes   = Var1_n_genes,
    mapped_icd10_range = icd10_category
  ) %>% 
  mutate(
    grSignature = case_when(
      str_detect(grSignature, "BloodCellsUp") ~ "bloodUp",
      str_detect(grSignature, "BloodCellsDown") ~ "bloodDown",
      str_detect(grSignature, "LungCellsUp") ~ "lungUp",
      str_detect(grSignature, "LungCellsDown") ~ "lungDown",
      str_detect(grSignature, "NeuralCellsUp") ~ "neuralUp",
      str_detect(grSignature, "NeuralCellsDown") ~ "neuralDown",
      str_detect(grSignature, "global_GR_genes_globalUp") ~ "systemicUp",
      str_detect(grSignature, "global_GR_genes_globalDown") ~ "systemicDown",
      TRUE ~ grSignature
    )
  ) %>% 
  mutate(
    # ---- regulation ----
    regulation = case_when(
      str_detect(grSignature, "Up") ~ "up",
      str_detect(grSignature, "Down") ~ "down",
      TRUE ~ NA_character_
    ),
    
    # ---- tissue ----
    tissue = case_when(
      str_detect(grSignature, "blood") ~ "blood",
      str_detect(grSignature, "lung") ~ "lung",
      str_detect(grSignature, "neural") ~ "neural",
      str_detect(grSignature, "systemic") ~ "systemic",
      TRUE ~ NA_character_
    ),
    
    # ---- signature type ----
    signatureType = case_when(
      tissue == "systemic" ~ "systemic",
      tissue %in% c("blood","lung","neural") ~ "tissue",
      TRUE ~ NA_character_
    )
  ) %>% 
  select(
    mapped_icd10_range,
    phenotype,
    phenotype_nGenes,
    grSignature,
    tissue,
    signatureType,
    regulation,
    grSignature_nGenes,
    expected_overlap,
    observed_overlap,
    odds_ratio,
    chi2,
    p_value,
    combine_score,
    overlap_genes,
    matched,
    n_unique_patterns
  ) %>% 
  select(-matched, -n_unique_patterns)

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

# ---- desired order of signatures within sheet ----
grSignature_order <- c(
  "systemicUp","systemicDown",
  "neuralUp","neuralDown",
  "bloodUp","bloodDown",
  "lungUp","lungDown"
)


phenotype_mapper <- c(
  # F3x — mood (affective) disorders
  "DisGeNET_Bipolar_Disorder" = "bipolar disorder",
  "DisGeNET_Depressive_disorder" = "depressive disorder",
  "DisGeNET_Major_Depressive_Disorder" = "major depressive disorder",
  "DisGeNET_Depression,_Bipolar" = "bipolar disorder - depressive episode",
  "DisGeNET_Bipolar_I_disorder,_most_recent_episode_manic_(disorder)" = "bipolar disorder - manic episode",
  "DisGeNET_Mixed_bipolar_I_disorder" = "mixed bipolar disorder",
  
  # F2x / psychotic disorders
  "DisGeNET_Psychotic_Disorders" = "psychotic disorders",
  "DisGeNET_Schizophrenia" = "schizophrenia",
  "DisGeNET_Nonorganic_psychosis" = "nonorganic psychosis",
  "genebass_Date_F23_first_reported_(acute_and_transient_psychotic_disorders)" = "acute psychotic disorders",
  "GWASCatalog_schizophrenia,intelligence,self reported educational attainment" = "schizophrenia & cognitive traits",
  "DisGeNET_Paranoia" = "paranoia",
  
  # F0x / organic, dementia, delirium
  "genebass_Date_F05_first_reported_(delirium,_not_induced_by_alcohol_and_other_psychoactive_substances)" = "non-substance-induced delirium",
  "DisGeNET_Delirium,_Dementia,_Amnestic,_Cognitive_Disorders" = "delirium, dementia & cognition",
  "GWASCatalog_dementia" = "dementia",
  "genebass_Date_F06_first_reported_(other_mental_disorders_due_to_brain_damage_and_dysfunction_and_to_physical_disease)" = "brain damage-related mental disorder",
  "genebass_Date_F01_first_reported_(vascular_dementia)" = "vascular dementia",
  "genebass_Date_F00_first_reported_(dementia_in_alzheimer's_disease)" = "dementia in Alzheimer's disease",
  
  # F1x / substance use
  "DisGeNET_Alcoholic_Intoxication,_Chronic" = "chronic alcohol intoxication",
  "DisGeNET_Alcoholic_Intoxication" = "alcohol intoxication",
  "genebass_Date_F17_first_reported_(mental_and_behavioural_disorders_due_to_use_of_tobacco)" = "tobacco use disorder",
  "genebass_Date_F11_first_reported_(mental_and_behavioural_disorders_due_to_use_of_opioids)" = "opioid use disorder",
  "genebass_Frequency_of_inability_to_cease_drinking_in_last_year" = "impaired control over drinking",
  "genebass_Age_when_last_took_cannabis" = "age at last cannabis use",
  
  # F4x / anxiety, neurotic, dissociative
  "genebass_Date_F48_first_reported_(other_neurotic_disorders)" = "other neurotic disorders",
  "genebass_Date_F41_first_reported_(other_anxiety_disorders)" = "other anxiety disorder",
  "genebass_Recent_feelings_of_foreboding" = "anxious foreboding",
  "genebass_Date_F44_first_reported_(dissociative_[conversion]_disorders)" = "dissociative disorders",
  "genebass_Recent_easy_annoyance_or_irritability" = "easy annoyance or irritability",
  "genebass_Date_F40_first_reported_(phobic_anxiety_disorders)" = "phobic anxiety disorders",
  
  # F8x / developmental disorders
  "DisGeNET_Learning_Disorders" = "learning disorders",
  "DisGeNET_Learning_Disabilities" = "learing disabilities",
  "DisGeNET_Autism_Spectrum_Disorders" = "autism spectrum disorders",
  "genebass_Date_F80_first_reported_(specific_developmental_disorders_of_speech_and_language)" = "speech & language disorders",
  "genebass_Date_F81_first_reported_(specific_developmental_disorders_of_scholastic_skills)" = "scholastic skill disorders",
  "DisGeNET_Developmental_Disabilities" = "developmental disorders"
)

# ---- DATA PIPELINE ----
df_final <- df_in %>%
  mutate(
    mapped_icd10_range = icd10_relabel(mapped_icd10_range)
  ) %>%
  filter(!is.na(mapped_icd10_range)) %>%
  mutate(
    grSignature = factor(grSignature, levels = grSignature_order),
    mapped_icd10_range = factor(mapped_icd10_range, levels = icd10_sheet_order)
  ) %>% 
  mutate(
    source    = str_extract(phenotype, "^[^_]+"),
    phenotype = str_remove(phenotype, "^[^_]+_F[0-9]x_")
  ) %>% 
  mutate(phenotype_v2 = paste0(source, "_", phenotype)) %>% 
  mutate(
    phenotype_label = dplyr::coalesce(unname(phenotype_mapper[phenotype_v2]), "-")
  ) %>% 
  select(
    mapped_icd10_range,
    phenotype,
    phenotype_label,
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
  



# ---- OUTPUT PATH ----
out_xlsx <- "/home/mateusz/projects/ifpan-michkor-gr/results_v2/overlap/overlap_grTissues_GWASCatalogDisGeNETgenebass/tables/supplementFile_overlap_grSystemicTissues_v2_06.03.2026.xlsx"

# ---- workbook + styles ----
wb <- createWorkbook()

# styl podświetlenia + border
highlight_style <- createStyle(
  fgFill = "#C9CBA3",
  border = "TopBottomLeftRight",
  borderColour = "#D9D9D9"
)

# styl zwykłych komórek z borderem
body_border_style <- createStyle(
  border = "TopBottomLeftRight",
  borderColour = "#D9D9D9"
)

# styl nagłówka z borderem
header_style <- createStyle(
  textDecoration = "bold",
  border = "TopBottomLeftRight",
  borderColour = "#D9D9D9"
)

# number formats
fmt_3dec <- createStyle(
  numFmt = "0.000",
  border = "TopBottomLeftRight",
  borderColour = "#D9D9D9"
)

fmt_sci3 <- createStyle(
  numFmt = "0.000E+00",
  border = "TopBottomLeftRight",
  borderColour = "#D9D9D9"
)

for (lvl in icd10_sheet_order) {
  
  df_sheet <- df_final %>%
    filter(mapped_icd10_range == lvl) %>%
    arrange(grSignature, desc(combine_score)) %>%
    mutate(grSignature = as.character(grSignature))
  
  if (nrow(df_sheet) == 0) next
  
  addWorksheet(wb, sheetName = lvl)
  writeData(wb, sheet = lvl, x = df_sheet, withFilter = TRUE)
  
  freezePane(wb, sheet = lvl, firstRow = TRUE)
  setColWidths(wb, sheet = lvl, cols = 1:ncol(df_sheet), widths = "auto")
  
  # ---- header style ----
  addStyle(
    wb, sheet = lvl, style = header_style,
    rows = 1, cols = 1:ncol(df_sheet),
    gridExpand = TRUE, stack = TRUE
  )
  
  # ---- base borders for all data cells ----
  addStyle(
    wb, sheet = lvl, style = body_border_style,
    rows = 2:(nrow(df_sheet) + 1),
    cols = 1:ncol(df_sheet),
    gridExpand = TRUE, stack = TRUE
  )
  
  # ---- highlight significant rows ----
  highlight_rows <- which(df_sheet$p_value < 0.05 & df_sheet$observed_overlap >= 3)
  if (length(highlight_rows) > 0) {
    addStyle(
      wb, sheet = lvl, style = highlight_style,
      rows = highlight_rows + 1,  # +1 for header
      cols = 1:ncol(df_sheet),
      gridExpand = TRUE, stack = TRUE
    )
  }
  
  # ---- apply number formats ----
  data_rows <- 2:(nrow(df_sheet) + 1)
  
  col_expected <- match("expected_overlap", names(df_sheet))
  col_cs       <- match("combine_score", names(df_sheet))
  
  col_chi2     <- match("chi2", names(df_sheet))
  col_p        <- match("p_value", names(df_sheet))
  col_or       <- match("odds_ratio", names(df_sheet))
  
  # expected_overlap & combine_score -> 3 decimals
  cols_3dec <- c(col_expected, col_cs)
  cols_3dec <- cols_3dec[!is.na(cols_3dec)]
  
  if (length(cols_3dec) > 0) {
    addStyle(
      wb, lvl, fmt_3dec,
      rows = data_rows, cols = cols_3dec,
      gridExpand = TRUE, stack = TRUE
    )
  }
  
  # chi2, p_value, odds_ratio -> scientific, 3 decimals
  cols_sci <- c(col_chi2, col_p, col_or)
  cols_sci <- cols_sci[!is.na(cols_sci)]
  
  if (length(cols_sci) > 0) {
    addStyle(
      wb, lvl, fmt_sci3,
      rows = data_rows, cols = cols_sci,
      gridExpand = TRUE, stack = TRUE
    )
  }
}

saveWorkbook(wb, out_xlsx, overwrite = TRUE)
message("Saved: ", out_xlsx)
# rm(df_in, df_final, out_xlsx, wb, df_sheet)



summary_df <- df_final %>% 
  mutate(
    signatureType = dplyr::case_when(
      grepl("blood", grSignature) ~ "blood",
      grepl("neural", grSignature) ~ "neural",
      grepl("lung", grSignature) ~ "lung", 
      grepl("systemic", grSignature) ~ "systemic"
    )
  ) %>% 
  select(-grSignature) %>% 
  group_by(mapped_icd10_range, signatureType) %>% 
  nest() %>% 
  mutate(category_phenotype_count = map_int(data, nrow)) %>%
  mutate(category_phenotype_count = map(data, ~ .x %>% .$phenotype %>% unique %>% length)) %>% unnest(category_phenotype_count) %>% 
  mutate(data_signif = map(data, ~ .x %>% filter(p_value < 0.05 & observed_overlap >= 3))) %>% 
  mutate(significant_association_count = map_int(data_signif, nrow)) %>% 
  mutate(significant_overlapping_gene_count = map_int(
    data_signif,
    ~ .x %>% .$overlap_genes %>% strsplit(",") %>% unlist() %>% unique() %>% length()
  )) %>% 
  mutate(significant_overlapping_genes = map_chr(
    data_signif,
    ~ .x %>% .$overlap_genes %>% strsplit(",") %>% unlist() %>% unique() %>% paste(collapse = ", ")
  )) %>% 
  mutate(combine_score = map_chr(data_signif, ~ paste(.x$combine_score, collapse = "|"))) %>% 
  mutate(signif_mean_combined_score = map_dbl(data_signif, ~ mean(.x$combine_score))) %>% 
  mutate(
    category_gene_count = case_when(
      mapped_icd10_range == "F30-F39" ~ 7558,
      mapped_icd10_range == "F20-F29" ~ 5555,
      mapped_icd10_range == "F00-F09" ~ 4474,
      mapped_icd10_range == "F10-F19" ~ 6157,
      mapped_icd10_range == "F40-F49" ~ 4973,
      mapped_icd10_range == "F80-F89" ~ 2834,
      mapped_icd10_range == "F50-F59" ~ 1099,
      mapped_icd10_range == "F60-F69" ~ 1539,
      mapped_icd10_range == "F70-F79" ~ 1381,
      mapped_icd10_range == "F90-F98" ~ 3295
    )
  ) %>% 
  mutate(
    category = case_when(
      mapped_icd10_range == "F30-F39" ~ "mood (affective) disorders",
      mapped_icd10_range == "F20-F29" ~ "schizophrenia, schizotypal and delusional disorders",
      mapped_icd10_range == "F00-F09" ~ "organic, including symptomatic, mental disorders",
      mapped_icd10_range == "F10-F19" ~ "mental & behavioural disorders due to psychoactive substance use",
      mapped_icd10_range == "F40-F49" ~ "neurotic, stress-related and somatoform disorders",
      mapped_icd10_range == "F80-F89" ~ "disorders of psychological development",
      mapped_icd10_range == "F50-F59" ~ "behavioural syndromes linked to physiological & somatic factors",
      mapped_icd10_range == "F60-F69" ~ "disorders of adult personality & behaviour",
      mapped_icd10_range == "F70-F79" ~ "mental retardation",
      mapped_icd10_range == "F90-F98" ~ "behavioural and emotional, childhood & adolescence onset"
    )
  ) %>% 
  select(-c(data, data_signif)) %>% 
  mutate(
    signatureType = factor(signatureType,
                           # levels = c("systemic", "neural", "blood", "lung"))
    levels = c("lung", "blood", "neural", "systemic"))
  ) %>% 
  select(
    category,
    mapped_icd10_range,
    signatureType,
    signif_mean_combined_score,
    significant_association_count,
    category_phenotype_count,
    significant_overlapping_gene_count,
    category_gene_count,
    significant_overlapping_genes,
    combine_score
  ) %>% 
  mutate(
    signif_mean_combined_score = ifelse(is.nan(signif_mean_combined_score), 0, signif_mean_combined_score)
  )

# ---- zapis do xlsx ----

file_path <- "/home/mateusz/projects/ifpan-michkor-gr/results_v2/overlap/overlap_grTissues_GWASCatalogDisGeNETgenebass/tables/summaryTable_grSystemicTissues_icd10F_v2_06.03.2026.xlsx"

wb <- createWorkbook()

addWorksheet(wb, "summary")

writeData(wb, "summary", summary_df)

# zamrożenie pierwszego wiersza
freezePane(wb, "summary", firstRow = TRUE)

saveWorkbook(wb, file_path, overwrite = TRUE)
  

df_final -> cleanOverlapResults_grSystemicTissues_ICD10F
# ##############################################################################
# ---- save to XLSX ----
# ##############################################################################
df_final %>% 
  mutate(
    signatureType = dplyr::case_when(
      grepl("blood", grSignature) ~ "blood",
      grepl("neural", grSignature) ~ "neural",
      grepl("lung", grSignature) ~ "lung", 
      grepl("systemic", grSignature) ~ "systemic"
    )
  ) %>% 
  select(-grSignature) %>% 
  select(-c(mapped_icd10_range, phenotype_label, regulation, expected_overlap, odds_ratio, chi2)) %>% 
  group_by(signatureType, source) %>% 
  nest() %>% 
  mutate(n_phenotypes = map(data, ~ .x %>% nrow)) %>% 
  unnest(n_phenotypes) %>% 
  mutate(n_signif_phenotypes = map(data, ~ .x %>% filter(p_value < 0.05) %>% filter(observed_overlap >= 3) %>% .$phenotype %>% unique %>% length)) %>% 
  unnest(n_signif_phenotypes) %>% 
  mutate(n_signif_genes = map(data, ~ .x %>% filter(p_value < 0.05) %>% filter(observed_overlap >= 3) %>% .$overlap_genes %>% strsplit(",") %>% unlist %>% unique %>% length)) %>% 
  unnest(n_signif_genes) %>% 
  mutate(
    n_all_genes = case_when(
      source == "DisGeNET" ~ 3395,
      source == "genebass" ~ 13313,
      source == "GWASCatalog" ~ 4114
    )
  ) %>% 
  mutate(prop_signif_genes = n_signif_genes/n_all_genes) %>% 
  mutate(prop_signif_phenotypes = n_signif_phenotypes/n_phenotypes) %>%
  mutate(n_ns_phenotypes = n_phenotypes - n_signif_phenotypes) %>% 
  mutate(n_ns_genes = n_all_genes - n_signif_genes) %>% 
  mutate(overlap_genes = map(data, ~ .x %>% filter(p_value < 0.05) %>% filter(observed_overlap >= 3) %>% .$overlap_genes %>% strsplit(",") %>% unlist %>% unique %>% paste(.,collapse = ", "))) %>% 
  select(c(signatureType, source, n_phenotypes, n_all_genes, n_signif_phenotypes, n_signif_genes, prop_signif_phenotypes, prop_signif_genes, overlap_genes)) %>% 
  unnest(overlap_genes) %>% 
  writexl::write_xlsx(
    "/home/mateusz/projects/ifpan-michkor-gr/results_v2/overlap/overlap_grTissues_GWASCatalogDisGeNETgenebass/tables/supplement_barplots_nAssociationsPerSource/summary_nPhenotypesGenes_perSignatureSource_v1_11.03.2026.xlsx"
  )


# ##############################################################################
df_final %>% 
  mutate(
    label = case_when(
      mapped_icd10_range == "F30-F39" ~ "F3x",
      mapped_icd10_range == "F20-F29" ~ "F2x",
      mapped_icd10_range == "F00-F09" ~ "F0x",
      mapped_icd10_range == "F10-F19" ~ "F1x",
      mapped_icd10_range == "F40-F49" ~ "F4x",
      mapped_icd10_range == "F80-F89" ~ "F8x",
      mapped_icd10_range == "F50-F59" ~ "F5x",
      mapped_icd10_range == "F60-F69" ~ "F6x",
      mapped_icd10_range == "F70-F79" ~ "F7x",
      mapped_icd10_range == "F90-F98" ~ "F9x"
    )
  ) %>% 
  mutate(
    label = paste0(source, "_", label, "_", phenotype)
  ) %>% 
  filter(source == "DisGeNET") %>% .$label %>% 
  unique() -> tmp

gene_list_disgenet[tmp] %>% unname() %>% unlist %>% unique %>% length()
gene_list_genebass[tmp] %>% unname %>% unlist %>% unique %>% length()
gene_list_GWASCatalog[tmp] %>% lapply(., length) %>% unname() %>% unlist %>% summary
