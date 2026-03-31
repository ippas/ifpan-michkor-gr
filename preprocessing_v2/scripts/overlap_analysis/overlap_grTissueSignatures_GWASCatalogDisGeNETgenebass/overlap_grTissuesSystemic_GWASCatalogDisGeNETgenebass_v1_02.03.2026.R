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


GWASCatalogDisGeNETgenebass_GrTissueSystemic_overlapChi2 <- run_full_overlap_analysis(
  gene_lists = c(flat_allGrSignatures_31.10.2025[sig_names], 
                 gene_list_disgenet,
                 gene_list_genebass,
                 gene_list_GWASCatalog
  ),
  total_genes = hgnc_symbols_vector_v110,
  cols_to_filter = c(sig_names),
  rows_to_filter = c(names(gene_list_disgenet),
                     names(gene_list_GWASCatalog),
                     names(gene_list_genebass)),
  plot_title_or = "",
  triangle_mode = "full",
  fdr_threshold = 1,
  data_type = "original_data",
  verbose = FALSE,
  palette_or = c("#c6d3e3", "white", "darkred"),
  text_contrast_range_or = c(-30, 4.9)
)

# adding score G
GWASCatalogDisGeNETgenebass_GrTissueSystemic_overlapChi2$processed$original_data$list$log10_number_overlap_matrix <- log10(GWASCatalogDisGeNETgenebass_GrTissueSystemic_overlapChi2$processed$original_data$list$number_overlap_matrix + 1)
GWASCatalogDisGeNETgenebass_GrTissueSystemic_overlapChi2$processed$original_data$list$log_p_value_matrix <- -log10(GWASCatalogDisGeNETgenebass_GrTissueSystemic_overlapChi2$processed$original_data$list$p_value_matrix) 


GWASCatalogDisGeNETgenebass_GrTissueSystemic_overlapChi2$processed$original_data$list$score_G <-
  GWASCatalogDisGeNETgenebass_GrTissueSystemic_overlapChi2$processed$original_data$list$log10_number_overlap_matrix *
  GWASCatalogDisGeNETgenebass_GrTissueSystemic_overlapChi2$processed$original_data$list$log_p_value_matrix 


GWASCatalogDisGeNETgenebass_GrTissueSystemic_overlapChi2$processed$original_data$df %>% head

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
  data_list = GWASCatalogDisGeNETgenebass_GrTissueSystemic_overlapChi2$processed,
  data_type = "original_data",
  
  # 🔺 triangle mode pełny
  triangle_mode = "full",
  
  
  color_key = "score_G",
  
  # 🔹 filtrowanie
  rows_to_filter = phenotypes,
  cols_to_filter = c(
    "global_GR_genes_globalUp5TissuesDerivedCells",
    "global_GR_genes_globalDown5TissuesDerivedCells",
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
    "global_GR_genes_globalUp5TissuesDerivedCells" = "systemicUp",
    "global_GR_genes_globalDown5TissuesDerivedCells" = "systemicDown",
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
  save_to_svg = "/home/mateusz/projects/ifpan-michkor-gr/results_v2/overlap/overlap_grTissues_GWASCatalogDisGeNETgenebass/figures/heatmap_NeuralBloodLungSystemic_top3GrSignatures_combineScore_v5_02.03.2026.svg",
  svg_width = 9.125,
  svg_height = 22,
  row_names_width = unit(90, "mm"),
  col_names_height = unit(50, "mm"),
  force_create_directory = TRUE
)


