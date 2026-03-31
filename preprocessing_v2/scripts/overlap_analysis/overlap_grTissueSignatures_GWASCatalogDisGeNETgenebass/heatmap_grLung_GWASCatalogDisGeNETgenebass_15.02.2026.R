GWASCatalogDisGeNETgenebass_GrTissue_overlapChi2$processed$original_data$df %>% 
  filter(grepl("Lung", Var2)) %>% 
  filter(p_value < 0.05) %>% 
  filter(gene_overlap_count >= 3) %>% 
  filter(!(Var1 %in% c("DisGeNET_NA", "genebass_NA", "GWASCatalog_NA",
                       "genebass_F9x_Date_F99_first_reported_(mental_disorder,_not_otherwise_specified)",
                       "GWASCatalog_F9x_F3x_F8x_F2x_attention deficit hyperactivity disorder,bipolar disorder,autism spectrum disorder,schizophrenia,major depressive disorder"))) %>% 
  filter(!grepl("neurodegenerative", Var1)) %>% 
  .$Var1 %>% unique -> phenotypes


GWASCatalogDisGeNETgenebass_GrTissue_overlapChi2$processed$original_data$df %>% 
  filter(grepl("Lung", Var2)) %>% 
  mutate(log10_pvalue = -log10(p_value)) %>% 
  mutate(log10_geneOverlap = log10(gene_overlap_count + 1)) %>% 
  mutate(combine_score = log10_pvalue * log10_geneOverlap) %>% 
  arrange(desc(combine_score)) %>% 
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
    icd10_category = map_chr(matched, ~ if(length(.x)==1) .x else NA_character_)
    # ---- zamiana nazw ----
  ) %>% 
  group_by(Var2, icd10_category) %>% 
  slice_max(order_by = combine_score, n = 4, with_ties = FALSE) %>%
  # filter(icd10_category %in% c("F0x", "F1x", "F2x", "F3x", "F4x", "F8x")) %>%
  ungroup() %>% 
  # filter(icd10_category ==  "F2x") %>% 
  # filter(p_value < 0.05) %>% 
  # filter(gene_overlap_count < 3) %>% 
  # filter(icd10_category == "F9x") %>% 
  # .$Var1
  group_by(icd10_category) %>% 
  # filter(combine_score > 0) %>% 
  filter(icd10_category == "F4x") %>% select(Var1, combine_score)
  select(Var1, icd10_category) %>%
  unique() %>% 
  # head %>% 
  arrange(icd10_category) %>% 
  .$Var1 -> phenotypes

phenotypes <- c(
  # F0x
  "DisGeNET_F0x_Other_specified_senile_psychotic_conditions",                                                 
  "genebass_F0x_Date_F00_first_reported_(dementia_in_alzheimer's_disease)",                                   
  "genebass_F0x_Date_F05_first_reported_(delirium,_not_induced_by_alcohol_and_other_psychoactive_substances)",
  
  "GWASCatalog_F0x_APOE carrier status,dementia",                                                            
  "GWASCatalog_F0x_cerebral amyloid angiopathy",                                                             
  "GWASCatalog_F0x_vascular dementia",    
  
  # F1x
  "genebass_F1x_Frequency_of_inability_to_cease_drinking_in_last_year",                                  
  "DisGeNET_F1x_Cannabis_Abuse",                                                                         
  "genebass_F1x_Ever_had_known_person_concerned_about,_or_recommend_reduction_of,_alcohol_consumption",  
  
  "DisGeNET_F1x_Marijuana_Abuse",                                                                        
  "DisGeNET_F1x_Alcoholic_Intoxication,_Chronic",                                                
  "DisGeNET_F1x_Cannabis_Dependence", 
  
  
  # F2x
  "genebass_F2x_Date_F25_first_reported_(schizoaffective_disorders)",               
  "genebass_F2x_Ever_believed_in_un-real_communications_or_signs",                  
  "genebass_F2x_Date_F22_first_reported_(persistent_delusional_disorders)",  
  
  "DisGeNET_F2x_Schizophrenia",                                                     
  "GWASCatalog_F2x_schizophrenia",                                                  
  "genebass_F2x_Date_F20_first_reported_(schizophrenia)", 
  
  # F3x
  "genebass_F3x_Recent_changes_in_speed/amount_of_moving_or_speaking",                     
  "GWASCatalog_F3x_I1x_mean arterial pressure,major depressive disorder",                  
  "genebass_F3x_Date_F39_first_reported_(unspecified_mood_[affective]_disorder)",         
  
  "GWASCatalog_F3x_major depressive disorder",                                             
  "DisGeNET_F3x_Major_Depressive_Disorder",                                                
  "DisGeNET_F3x_Unipolar_Depression",
  
  # F4x
  "genebass_F4x_Date_F45_first_reported_(somatoform_disorders)"  ,                                         
  "genebass_F4x_Recent_feelings_of_foreboding",                                           
  "genebass_F4x_Recent_inability_to_stop_or_control_worrying",
  
  "genebass_F4x_Date_F40_first_reported_(phobic_anxiety_disorders)",
  "GWASCatalog_F4x_anxiety disorder",                               
  "GWASCatalog_F4x_neurotic disorder",
  
  # F8x
  "DisGeNET_F8x_Neurodevelopmental_Disorders",                                                     
  "DisGeNET_F8x_Developmental_Coordination_Disorder",                                             
  "genebass_F8x_Date_F80_first_reported_(specific_developmental_disorders_of_speech_and_language)",
  
  "GWASCatalog_F8x_dyslexia",                                                                      
  "DisGeNET_F8x_Autism_Spectrum_Disorders",                                                        
  "DisGeNET_F8x_Learning_Disorders",   
  
  # F9x
  "GWASCatalog_F9x_attention deficit hyperactivity disorder",                      
  "DisGeNET_F9x_Gilles_de_la_Tourette_syndrome",                                   
  "DisGeNET_F9x_Stereotypic_Movement_Disorder", 
  
  # other
  "GWASCatalog_F5x_insomnia",                                                                                                        
  "genebass_F5x_Date_F52_first_reported_(sexual_dysfunction,_not_caused_by_organic_disorder_or_disease)" ,                           
  "DisGeNET_F5x_Erectile_dysfunction" ,                                                                                              
  "GWASCatalog_F5x_bulimia nervosa",                                                                                                 
  "genebass_F5x_Date_F50_first_reported_(eating_disorders)" ,                                                                        
  "genebass_F5x_Date_F53_first_reported_(mental_and_behavioural_disorders_associated_with_the_puerperium,_not_elsewhere_classified)",
  "GWASCatalog_F5x_eating disorder",
  
  "DisGeNET_F7x_Intellectual_Disability" ,                      
  "DisGeNET_F7x_Severe_intellectual_disability_(disorder)",     
  "genebass_F7x_Attempted_fluid_intelligence_(FI)_test.",       
  "DisGeNET_F7x_Non-specific_syndromic_intellectual_disability",
  "genebass_F7x_Fluid_intelligence_score",                      
  "DisGeNET_F7x_Coffin-Siris_syndrome",                         
  "DisGeNET_F7x_Mental_Retardation",
  "GWASCatalog_F6x_childhood gender nonconformity"
)

heatmap_overlap_log2OR_complex(
  data_list = GWASCatalogDisGeNETgenebass_GrTissue_overlapChi2$processed,
  data_type = "original_data",
  
  # 🔺 triangle mode pełny
  triangle_mode = "full",
  
  
  color_key = "score_G",
  
  # 🔹 filtrowanie
  rows_to_filter = phenotypes,
  cols_to_filter = c(
    "minusGlobalUpDown5TissuesDerivedCells_LungCellsUp",
    "minusGlobalUpDown5TissuesDerivedCells_LungCellsDown"
  ),
  row_order_original = phenotypes,
  
  # 🎨 skala kolorów
  color_scale_range = c(0, 10),
  text_contrast_range = c(-30, 5),
  palette = c("white", "#f0cfc7ff", "#bd6651ff","#992314ff", "darkred"),
  
  # 🔲 prostokąty istotności
  color_rects = c("#97C426", "#2F4603", "purple"),
  
  # 📊 klastrowanie
  cluster_rows = F,
  cluster_cols = F,
  show_dendrograms = TRUE,
  rect_lwd = 2.5,
  p_thresholds = c(0.05, 0.01, 0.001),
  col_mapper = c(
    "minusGlobalUpDown5TissuesDerivedCells_NeuralCellsUp" =  "neuralUp",
    "minusGlobalUpDown5TissuesDerivedCells_NeuralCellsDown" =  "neuralDown",
    "minusGlobalUpDown5TissuesDerivedCells_BloodCellsUp" = "bloodUp",
    "minusGlobalUpDown5TissuesDerivedCells_BloodCellsDown" = "bloodDown",
    "minusGlobalUpDown5TissuesDerivedCells_LungCellsUp" = "lungUp",
    "minusGlobalUpDown5TissuesDerivedCells_LungCellsDown" = "lungDown"
  ),
  row_dend_height = unit(20, "mm"),
  col_dend_height = unit(20, "mm"),
  tile_gap = 1,
  save_to_svg = "/home/mateusz/projects/ifpan-michkor-gr/results_v2/overlap/overlap_grTissues_GWASCatalogDisGeNETgenebass/figures/heatmap_grLung_GWASCatalogDisGeNETgenebass_combineScore_v3_23.02.2026.svg",
  svg_width = 6.625,
  svg_height = 22,
  row_names_width = unit(90, "mm"),
  col_names_height = unit(50, "mm"),
  force_create_directory = TRUE
)

# ##############################################################################
# ---- prepare supplement file ----
# ##############################################################################
GWASCatalogDisGeNETgenebass_GrTissue_overlapChi2$processed$original_data$df %>% 
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
                  "global_GR_genes_globalDown5TissuesDerivedCells" = "systemicDown", 
                  "minusGlobalUpDown5TissuesDerivedCells_BloodCellsUp" = "bloodUp",
                  "minusGlobalUpDown5TissuesDerivedCells_LungCellsUp" = "lungUp",
                  "minusGlobalUpDown5TissuesDerivedCells_NeuralCellsUp" = "neuralUp",
                  "minusGlobalUpDown5TissuesDerivedCells_BloodCellsDown" = "bloodDown",
                  "minusGlobalUpDown5TissuesDerivedCells_LungCellsDown" = "lungDown",
                  "minusGlobalUpDown5TissuesDerivedCells_NeuralCellsDown" = "neuralDown"
    )
  ) %>% 
  filter(Var2 %in% c("lungUp", "lungDown")) %>% 
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
                         levels = c("lungUp", "lungDown"))
  ) %>% 
  arrange(grSignature, desc(combine_score)) %>% 
  split(.$icd10_category) -> grLung_GWASCatalogDisGeNETgenebass_list


# ##############################################################################
# ---- save to file ----
# ##############################################################################
# -------------------------
# input
# -------------------------
df_list <- grLung_GWASCatalogDisGeNETgenebass_list

out_dir  <- "/home/mateusz/projects/ifpan-michkor-gr/results_v2/overlap/overlap_grTissues_GWASCatalogDisGeNETgenebass/tables/"
out_file <- file.path(out_dir, "overlap_grLung_GWASCatalogDisGeNETgenebassICD10F_v1_22.02.2026.xlsx")

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
grLung_GWASCatalogDisGeNETgenebass_list$F9x %>% 
  filter(p_value < 0.05) %>% 
  filter(observed_overlap >= 3) %>% 
  .$combine_score %>% mean

grLung_GWASCatalogDisGeNETgenebass_list$F9x %>% 
  filter(p_value < 0.05) %>% 
  filter(observed_overlap >= 3) %>% nrow


grLung_GWASCatalogDisGeNETgenebass_list$F9x %>% 
  filter(p_value < 0.05) %>% 
  filter(observed_overlap >= 3) %>% 
  .$overlap_genes %>% 
  strsplit(",") %>% unlist %>% unique %>% length()

grLung_GWASCatalogDisGeNETgenebass_list$F9x %>% 
  filter(p_value < 0.05) %>% 
  filter(observed_overlap >= 3) %>% 
  .$overlap_genes %>% 
  strsplit(",") %>% unlist %>% unique %>% 
  paste(collapse = ", ")

grLung_GWASCatalogDisGeNETgenebass_list$F9x %>% 
  filter(p_value < 0.05) %>% 
  filter(observed_overlap >= 3) %>% 
  .$combine_score %>% 
  paste(collapse = "|")

