GWASCatalogDisGeNETgenebass_GrTissue_overlapChi2$processed$original_data$df %>% 
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
      grepl("Up$", Var2)   ~ "Up",
      grepl("Down$", Var2) ~ "Down",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(regulation), !is.na(icd10_category)) %>%
  filter(icd10_category == "F9x") %>% filter(p_value < 0.05) %>% filter(gene_overlap_count > 8) %>% .$overlap_genes
  
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
  # ---- F3x ----
  "genebass_F3x_Recent_changes_in_speed/amount_of_moving_or_speaking",
  "genebass_F3x_Recent_trouble_concentrating_on_things",              
  "genebass_F3x_Recent_feelings_of_depression",                       
  "DisGeNET_F3x_Depressive_disorder",                                 
  "DisGeNET_F3x_Unipolar_Depression",                                 
  "DisGeNET_F3x_Major_Depressive_Disorder",
  # ---- F2x ----
  "genebass_F2x_Ever_believed_in_un-real_communications_or_signs",         
  "genebass_F2x_Date_F22_first_reported_(persistent_delusional_disorders)",
  "genebass_F2x_Date_F20_first_reported_(schizophrenia)",                  
  "DisGeNET_F2x_Schizophrenia",                                            
  "GWASCatalog_F2x_schizophrenia",                                         
  "GWASCatalog_F2x_E1x_schizophrenia,type 2 diabetes mellitus",   
  
  # ---- F1x ----
  "genebass_F1x_Date_F10_first_reported_(mental_and_behavioural_disorders_due_to_use_of_alcohol)",     
  "genebass_F1x_Frequency_of_inability_to_cease_drinking_in_last_year",                                
  "genebass_F1x_Ever_had_known_person_concerned_about,_or_recommend_reduction_of,_alcohol_consumption",
  "GWASCatalog_F1x_cocaine dependence",                                                                
  "GWASCatalog_F1x_cannabis dependence",                                                               
  "DisGeNET_F1x_Marijuana_Abuse",
  # ---- F0x ----
  "GWASCatalog_F0x_APOE carrier status,dementia",                                          
  "GWASCatalog_F0x_dementia",                                                              
  "DisGeNET_F0x_Delirium,_Dementia,_Amnestic,_Cognitive_Disorders",                        
  "GWASCatalog_F0x_vascular dementia",                                                     
  "DisGeNET_F0x_Other_specified_senile_psychotic_conditions",                              
  "genebass_F0x_Date_F02_first_reported_(dementia_in_other_diseases_classified_elsewhere)",
  # ---- F4x ----
  "genebass_F4x_Recent_feelings_of_foreboding",                                
  "genebass_F4x_Recent_easy_annoyance_or_irritability",                        
  "genebass_F4x_Date_F42_first_reported_(obsessive-compulsive_disorder)",      
  "GWASCatalog_F4x_hoarding disorder",                                         
  "genebass_F4x_Date_F44_first_reported_(dissociative_[conversion]_disorders)",
  "genebass_F4x_Date_F40_first_reported_(phobic_anxiety_disorders)",  
  # ---- F5x ----
  "genebass_F5x_Date_F52_first_reported_(sexual_dysfunction,_not_caused_by_organic_disorder_or_disease)",
  "GWASCatalog_F5x_eating disorder",                                                                     
  "genebass_F5x_Date_F50_first_reported_(eating_disorders)",                                             
  "GWASCatalog_F5x_insomnia",                                                                            
  "DisGeNET_F5x_Erectile_dysfunction",                                                                   
  "genebass_F5x_Date_F51_first_reported_(nonorganic_sleep_disorders)", 
  # ---- F8x ----
  "genebass_F8x_Date_F80_first_reported_(specific_developmental_disorders_of_speech_and_language)",
  "DisGeNET_F8x_Autistic_Disorder",                                                                
  "DisGeNET_F8x_Neurodevelopmental_Disorders",                                                     
  "GWASCatalog_F8x_dyslexia",                                                                      
  "DisGeNET_F8x_Autism_Spectrum_Disorders",                                                        
  "DisGeNET_F8x_Learning_Disorders",     
  # ---- F9x ----
  "GWASCatalog_F9x_conduct disorder",                                              
  "GWASCatalog_F9x_Tourette syndrome",                                             
  "GWASCatalog_F9x_attention deficit hyperactivity disorder,conduct disorder",     
  "GWASCatalog_F9x_attention deficit hyperactivity disorder",                      
  "GWASCatalog_F9x_attention deficit hyperactivity disorder,risk-taking behaviour",
  "DisGeNET_F9x_Gilles_de_la_Tourette_syndrome",
  # ---- F9x ----
  "GWASCatalog_F9x_conduct disorder",                                              
  "GWASCatalog_F9x_Tourette syndrome",                                             
  "GWASCatalog_F9x_attention deficit hyperactivity disorder,conduct disorder",     
  "GWASCatalog_F9x_attention deficit hyperactivity disorder",                      
  "GWASCatalog_F9x_attention deficit hyperactivity disorder,risk-taking behaviour",
  "DisGeNET_F9x_Gilles_de_la_Tourette_syndrome"  
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
  save_to_svg = "/home/mateusz/projects/ifpan-michkor-gr/results_v2/overlap/overlap_grTissues_GWASCatalogDisGeNETgenebass/figures/heatmap_top3NeuralBloodLung_GWASCatalogDisGeNETgenebass_combineScore_v3_27.02.2026.svg",
  svg_width = 8.225,
  svg_height = 22,
  row_names_width = unit(90, "mm"),
  col_names_height = unit(50, "mm"),
  force_create_directory = TRUE
)

