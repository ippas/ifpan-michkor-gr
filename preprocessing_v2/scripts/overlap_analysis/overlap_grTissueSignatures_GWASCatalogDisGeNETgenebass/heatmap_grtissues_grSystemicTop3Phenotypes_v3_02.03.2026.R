GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$df %>% 
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
    icd10_category = map_chr(matched, ~ .x[[1]]),
    regulation = case_when(
      grepl("global_GR_genes_globalUp5TissuesDerivedCells$", Var2)   ~ "Up",
      grepl("global_GR_genes_globalDown5TissuesDerivedCells$", Var2) ~ "Down",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(regulation), !is.na(icd10_category)) %>%
  filter(icd10_category == "F9x") %>%
  
  # 1) SUMY per (fenotyp, regulacja)
  group_by(icd10_category, Var1, regulation) %>%
  summarise(sum_cs = sum(combine_score, na.rm = TRUE), .groups = "drop") %>%
  
  # 2) Up/Down w kolumnach
  pivot_wider(
    names_from = regulation,
    values_from = sum_cs,
    values_fill = 0
  ) %>%
  
  # 3) winner (unikat Var1) + tiebreak (remis -> Up, żeby nie gubić)
  mutate(
    winner = if_else(Up >= Down, "Up", "Down"),
    winner_score = pmax(Up, Down)
  ) %>%
  
  # 4) wybór TOP5 per Up/Down po winner, a jeśli brakuje -> dopełnij najlepszymi po tej regulacji
  group_by(icd10_category) %>%
  group_modify(~{
    df <- .x
    
    pick_side <- function(side, n = 5) {
      # najpierw winnerzy
      w <- df %>%
        filter(winner == side) %>%
        arrange(desc(winner_score)) %>%
        distinct(Var1, .keep_all = TRUE) %>%
        slice_head(n = n)
      
      # jeśli mało, dopełnij przegranymi o najwyższym sum_cs dla tej strony (bez duplikatów)
      if (nrow(w) < n) {
        need <- n - nrow(w)
        extra <- df %>%
          filter(!(Var1 %in% w$Var1)) %>%
          arrange(desc(.data[[side]])) %>%   # sortuj po Up albo Down
          distinct(Var1, .keep_all = TRUE) %>%
          slice_head(n = need)
        
        w <- bind_rows(w, extra)
      }
      
      w %>%
        mutate(winner = side) %>%
        mutate(rank_within = row_number())
    }
    
    bind_rows(
      pick_side("Up", 3),
      pick_side("Down", 3)
    )
  }) %>%
  ungroup() %>%
  mutate(winner = factor(winner, levels = c("Up", "Down"))) %>%
  arrange(winner, rank_within) %>% 
  pull(Var1) -> phenotypes

phenotypes



phenotypes <- c(
  # F0x
  "genebass_F0x_Date_F05_first_reported_(delirium,_not_induced_by_alcohol_and_other_psychoactive_substances)",
  "DisGeNET_F0x_Delirium,_Dementia,_Amnestic,_Cognitive_Disorders", 
  "GWASCatalog_F0x_dementia", 
  
  "genebass_F0x_Date_F06_first_reported_(other_mental_disorders_due_to_brain_damage_and_dysfunction_and_to_physical_disease)",
  "genebass_F0x_Date_F01_first_reported_(vascular_dementia)", 
  "genebass_F0x_Date_F00_first_reported_(dementia_in_alzheimer's_disease)",
  
  # F1x
  "DisGeNET_F1x_Alcoholic_Intoxication,_Chronic",                                                                             
  "DisGeNET_F1x_Alcoholic_Intoxication", 
  "genebass_F1x_Date_F17_first_reported_(mental_and_behavioural_disorders_due_to_use_of_tobacco)",                            
  
  "genebass_F1x_Frequency_of_inability_to_cease_drinking_in_last_year",                                                       
  "genebass_F1x_Date_F11_first_reported_(mental_and_behavioural_disorders_due_to_use_of_opioids)",                            
  "genebass_F1x_Age_when_last_took_cannabis",
  
  # F2x
  "DisGeNET_F2x_Psychotic_Disorders",
  "DisGeNET_F2x_Schizophrenia",
  "DisGeNET_F2x_Nonorganic_psychosis", 
  
  "genebass_F2x_Date_F23_first_reported_(acute_and_transient_psychotic_disorders)",
  "GWASCatalog_F2x_schizophrenia,intelligence,self reported educational attainment",    
  "DisGeNET_F2x_Paranoia",
  # "genebass_F2x_Date_F22_first_reported_(persistent_delusional_disorders)",  
  
  # F3x
  "DisGeNET_F3x_Bipolar_Disorder",                                                                                            
  "DisGeNET_F3x_Depressive_disorder",
  "DisGeNET_F3x_Major_Depressive_Disorder",
  
  "DisGeNET_F3x_Mixed_bipolar_I_disorder",                                                                                    
  "DisGeNET_F3x_Bipolar_I_disorder,_most_recent_episode_manic_(disorder)",                                                    
  "DisGeNET_F3x_Depression,_Bipolar",                                                                                         

  # F4x
  "genebass_F4x_Date_F48_first_reported_(other_neurotic_disorders)",                                                         
  "genebass_F4x_Date_F41_first_reported_(other_anxiety_disorders)",                                                           
  "genebass_F4x_Date_F44_first_reported_(dissociative_[conversion]_disorders)",
  
  "genebass_F4x_Recent_feelings_of_foreboding",                                                                               
  "genebass_F4x_Date_F40_first_reported_(phobic_anxiety_disorders)" ,                                                         
  "genebass_F4x_Recent_easy_annoyance_or_irritability" ,         

  # F8x
  "DisGeNET_F8x_Learning_Disorders",                                                                                          
  "DisGeNET_F8x_Learning_Disabilities",                                                                                      
  "DisGeNET_F8x_Autism_Spectrum_Disorders",
  
  "genebass_F8x_Date_F81_first_reported_(specific_developmental_disorders_of_scholastic_skills)",                             
  "DisGeNET_F8x_Developmental_Disabilities",                                                                                  
  "genebass_F8x_Date_F80_first_reported_(specific_developmental_disorders_of_speech_and_language)",
  
  # F5x
  "genebass_F5x_Date_F51_first_reported_(nonorganic_sleep_disorders)",                                                               
  "genebass_F5x_Date_F53_first_reported_(mental_and_behavioural_disorders_associated_with_the_puerperium,_not_elsewhere_classified)",
  "genebass_F5x_Date_F50_first_reported_(eating_disorders)",                                                                         
  "DisGeNET_F5x_Erectile_dysfunction",                                                                                               
  "genebass_F5x_Date_F52_first_reported_(sexual_dysfunction,_not_caused_by_organic_disorder_or_disease)",   
  "DisGeNET_F5x_Pediatric_failure_to_thrive",
  ############################3
  # F9x
  "DisGeNET_F9x_Gilles_de_la_Tourette_syndrome",
  "GWASCatalog_F9x_conduct disorder",
  "GWASCatalog_F9x_attention deficit hyperactivity disorder,conduct disorder",
  "GWASCatalog_F9x_attention deficit hyperactivity disorder",                      
  "DisGeNET_F9x_Stereotypic_Movement_Disorder",                                    
  "GWASCatalog_F9x_attention deficit hyperactivity disorder,risk-taking behaviour",
  
  # other
  "genebass_F5x_Date_F51_first_reported_(nonorganic_sleep_disorders)",                                                               
  "GWASCatalog_F9x_conduct disorder",
  "genebass_F5x_Date_F50_first_reported_(eating_disorders)",                                                                         
  "DisGeNET_F5x_Erectile_dysfunction",                                                                                               
  "genebass_F5x_Date_F52_first_reported_(sexual_dysfunction,_not_caused_by_organic_disorder_or_disease)",   
  "DisGeNET_F5x_Pediatric_failure_to_thrive"
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
    "minusGlobalUpDown5TissuesDerivedCells_NeuralCellsUp",
    "minusGlobalUpDown5TissuesDerivedCells_BloodCellsUp",
    "minusGlobalUpDown5TissuesDerivedCells_LungCellsUp",
    "minusGlobalUpDown5TissuesDerivedCells_NeuralCellsDown",
    "minusGlobalUpDown5TissuesDerivedCells_BloodCellsDown",
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
  save_to_svg = "/home/mateusz/projects/ifpan-michkor-gr/results_v2/overlap/overlap_grTissues_GWASCatalogDisGeNETgenebass/figures/heatmap_NeuralBloodLung_top3GrSignatures_combineScore_v4_02.03.2026.svg",
  svg_width = 8.225,
  svg_height = 22,
  row_names_width = unit(90, "mm"),
  col_names_height = unit(50, "mm"),
  force_create_directory = TRUE
)

