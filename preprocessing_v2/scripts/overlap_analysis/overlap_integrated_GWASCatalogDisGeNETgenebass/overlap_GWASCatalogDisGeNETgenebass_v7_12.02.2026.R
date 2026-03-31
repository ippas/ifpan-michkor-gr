# ##############################################################################
# ---- uses data ----
# ##############################################################################

# neurocognitive–neurodegenerative spectrum (F01–F09, G10, G30–G32, G35 – G37)
# substance use (F10 – F19), 
# psychotic spectrum (F20 – F29), 
# mood disorders (F30 – F39), 
# anxiety spectrum (F40 – F49), 
# physiological-behavioral syndromes (F50 – F59),
# personality-behavioral spectrum (F60 – F69),
# intellectual disabilities (F70 – F79), 
# neurodevelopmental spectrum (F80 – F89), 
# behavioral-emotional spectrum (F90 – F99)

sig_names <- c(
  # "minusGlobalUpDown5TissuesDerivedCells_BloodCellsUp",
  # "minusGlobalUpDown5TissuesDerivedCells_LungCellsUp",
  # "minusGlobalUpDown5TissuesDerivedCells_NeuralCellsUp",
  # "minusGlobalUpDown5TissuesDerivedCells_BloodCellsDown",
  # "minusGlobalUpDown5TissuesDerivedCells_LungCellsDown",
  # "minusGlobalUpDown5TissuesDerivedCells_NeuralCellsDown",
  "global_GR_genes_globalUp5TissuesDerivedCells",
  "global_GR_genes_globalDown5TissuesDerivedCells"
)


# mapping phenotypes to ICD10
mapperPhenotyeps2icd10 <- read_excel_sheets("data/mapperPhenotypesICD10/manualMapping_phenotype2ICD10categories_05.02.2026.xlsx") %>% 
  set_names(c( "GWASCatalog", "DisGeNET_score0", "genebass", "DisGeNET_score0.5"))

# ---- mapper GWAS Catalog phenotypes ----
mapperPhenotyeps2icd10$GWASCatalog %>% 
  mutate(
    mapping2ICD10_v2 = mapping2ICD10_v2 %>%
      str_replace_all("[ \\+\\(\\);]", "_") %>%  # spacja, +, (, ), ;
      str_replace_all("_+", "_") %>%             # wiele _ → jedno _
      str_replace_all("^_|_$", "")               # usuń _ na początku/końcu
  ) %>% 
  mutate(names = paste0(mapping2ICD10_v2, "_", phenotype)) %>% 
  filter(grepl("F|neurodegenerative", mapping2ICD10_v2)) %>% 
  select(c(phenotype, names)) -> mapperPhenotypesGWASCatalog

# ---- mapper genebass phenotypes ----
mapperPhenotyeps2icd10$genebass %>% 
  mutate(phenotype = str_replace_all(phenotype, " ", "_")) %>% 
  # mutate(names = paste0(subcategory, "_", phenotype)) %>% 
  # mutate(names = str_remove(names, "NA_")) %>% 
  mutate(names = phenotype) %>% 
  mutate(names = paste0(mapping2ICD10_v2, "_", names)) %>% 
  filter(grepl("F|neurodegenerative", mapping2ICD10_v2)) %>% 
  select(c(phenotype, names)) -> mapperPhenotypesGenebass

# ---- mapper DisGeNET phenotypes ----
mapperPhenotyeps2icd10$DisGeNET_score0 %>% 
  mutate(phenotype = str_replace_all(phenotype, " ", "_")) %>%
  mutate(
    mapping2ICD10_v2 = mapping2ICD10_v2 %>%
      str_replace_all("[ \\+\\(\\);]", "_") %>%  # spacja, +, (, ), ;
      str_replace_all("_+", "_") %>%             # wiele _ → jedno _
      str_replace_all("^_|_$", "")               # usuń _ na początku/końcu
  ) %>% 
  filter(grepl("F|neurodegenerative", mapping2ICD10_v2)) %>% 
  mutate(names = paste0(mapping2ICD10_v2, "_", phenotype)) %>% 
  select(c(phenotype, names)) -> mapperPhenotypesDisGeNET


# ##############################################################################
# ---- GWAS Catalog ----
# ##############################################################################
efo <- "EFO_0000677"
url_assoc <- sprintf(
  "https://www.ebi.ac.uk/gwas/api/v2/efotraits/%s/associations/download?includeBgTraits=false&includeChildTraits=true",
  efo
)

moodDisorders_GWASCatalog <- read_tsv(url_assoc, show_col_types = FALSE, progress = FALSE)

moodDisorders_GWASCatalog %>%
  select(mappedGenes, efoTraits, pValue) %>% 
  separate_rows(mappedGenes, sep = ",") %>% 
  filter(mappedGenes %in% hgnc_symbols_vector_v110) %>% 
  left_join(., mapperPhenotypesGWASCatalog, by = c("efoTraits" = "phenotype")) %>%
  # mutate(names = efoTraits) %>% 
  group_by(names) %>% 
  nest %>% 
  mutate(n_genes = map(data, ~ .x %>% .$mappedGenes %>% unique %>%  length)) %>% 
  unnest(n_genes) %>%
  filter(n_genes >= 10) %>% 
  unnest() %>% 
  select(c(names, mappedGenes)) %>% 
  mutate(names = paste0("GWASCatalog_", names)) %>% 
  { split(.$mappedGenes, .$names) } %>% 
  lapply(unique) -> gene_list_GWASCatalog

# ##############################################################################
# ---- genebass ----
# ##############################################################################
input_dir <- "data/genebass_v2/all_categories_SKAT/"

files <- c(
  list.files(
    path = input_dir,
    pattern = "SKAT_Online_follow-up_.*Mental_health.*\\.tsv\\.bgz$",
    full.names = TRUE
  ),
  file.path(input_dir, "genebass_SKAT_Health-related_outcomes_-_First_occurrences_-_Mental_and_behavioural_disorders.tsv.bgz"),
  file.path(input_dir, "genebass_SKAT_Health-related_outcomes_-_First_occurrences_-_Nervous_system_disorders.tsv.bgz"),
  file.path(input_dir, "genebass_SKAT_Health-related_outcomes_-_First_occurrences_-_Mental_and_behavioural_disorders.tsv.bgz"),
  file.path(input_dir, "genebass_SKAT_Health-related_outcomes_-_First_occurrences_-_Congenital_disruptions_and_chromosomal_abnormalities.tsv.bgz"),
  file.path(input_dir, "genebass_SKAT_UK_Biobank_Assessment_Centre_-_Cognitive_function_-_Fluid_intelligence___reasoning.tsv.bgz")
)


genebass_online_mental <- files %>%
  map_df(
    ~ read_tsv(.x, show_col_types = FALSE) %>%
      select(-c(pvalue_test, pvalue_threshold, heritability)) %>%
      filter(annotation == "pLoF", pvalue < 0.05) %>%
      mutate(
        source_file = basename(.x),
        description_format = str_replace_all(description, " ", "_"),
        source_trait = str_replace(source_file, "genebass_burden__-_", "") %>%
          str_replace(".tsv.bgz", "")
      )
  )


gene_list_genebass <- genebass_online_mental %>% 
  filter(pvalue < 0.01) %>% 
  filter(gene_symbol %in% hgnc_symbols_vector_v110) %>% 
  select(c(gene_symbol, description_format)) %>% 
  unique() %>% 
  group_by(description_format) %>% 
  nest %>% 
  mutate(n_genes = map(data, ~ .x %>% nrow)) %>% 
  unnest(n_genes) %>% 
  filter(n_genes >= 10) %>%
  select(-n_genes) %>% 
  unnest(data) %>% 
  ungroup %>% 
  distinct() %>%
  left_join(., mapperPhenotypesGenebass, by = c("description_format" = "phenotype")) %>% 
  select(-description_format) %>% 
  mutate(names = paste0("genebass_", names)) %>% 
  group_by(names) %>%
  summarise(genes = list(unique(gene_symbol)), .groups = "drop") %>%
  deframe()

# ##############################################################################
# ---- DisGeNET ----
# ##############################################################################

disgenet_mentalDisorders <- readRDS("data/databases/disgenet/disgenet_mentalDisordersF03.rds")

filter_min_vector_length <- function(x, min_len = 10) {
  Filter(function(v) length(v) >= min_len, x)
}


gene_list_disgenet <- disgenet_mentalDisorders$geneLists_scoreMin0 %>%
  filter_min_vector_length(min_len = 10)

# gene_list_disgenet <- disgenet_mentalDisorders$geneLists_scoreMin0.5 %>%
#   filter_min_vector_length(min_len = 10)

gene_list_disgenet %>% stack %>% 
  set_colnames(c("gene_symbol", "phenotype")) %>% 
  left_join(., mapperPhenotypesDisGeNET, by = "phenotype") %>% 
  select(names, gene_symbol) %>% 
  group_by(names) %>%
  mutate(names = paste0("DisGeNET_", names)) %>% 
  summarise(gene_symbol = list(gene_symbol), .groups = "drop") %>%
  deframe() -> gene_list_disgenet

# ##############################################################################
# ---- chi2 analysis ----
# ##############################################################################
GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2 <- run_full_overlap_analysis(
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

# ##############################################################################
# ---- selection phenotypes ----
# ##############################################################################


GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$df %>% 
  filter(grepl("F0x", Var1)) %>% 
  # filter(gene_overlap_count >= 3) %>% 
  mutate(
    source = case_when(
      Var1 %in% names(gene_list_disgenet) ~ "DisGeNET",
      Var1 %in% names(gene_list_genebass) ~ "genebass",
      Var1 %in% names(gene_list_GWASCatalog) ~ "GWAS_Catalog",
      TRUE        ~ "default_value"   # fallback / else
    )
  ) %>% 
  # mutate(Var1 = paste0(source, "_", Var1)) %>% 
  group_by(source) %>% 
  slice_min(order_by = p_value, n = 4, with_ties = FALSE) %>%
  arrange(source, p_value) %>% 
  # filter(grepl("GWAS", Var1)) %>%
  select(c(Var1, Var2, p_value, gene_overlap_count))
.$Var1 %>% unique
unique -> GWASCatalogDisGeNETgenebass_phenotypes


GWASCatalogDisGeNETgenebass_phenotypes

GWASCatalogDisGeNETgenebass_phenotypes <- c(
  # F0
  "DisGeNET_F0x_Delirium,_Dementia,_Amnestic,_Cognitive_Disorders",
  "DisGeNET_F0x_Delirium",
  "GWASCatalog_F0x_dementia",
  "GWASCatalog_F0x_vascular dementia",
  "genebass_F0x_Date_F05_first_reported_(delirium,_not_induced_by_alcohol_and_other_psychoactive_substances)",
  "genebass_F0x_Date_F06_first_reported_(other_mental_disorders_due_to_brain_damage_and_dysfunction_and_to_physical_disease)",
  # F1
  "DisGeNET_F1x_Alcoholic_Intoxication,_Chronic",
  "DisGeNET_F1x_Alcoholic_Intoxication",
  "GWASCatalog_F1x_nicotine use,generational effect measurement,illegal drug consumption,non-substance related disinhibited behaviour,alcohol drinking,alcohol dependence",
  "GWASCatalog_F1x_nicotine dependence",
  "genebass_F1x_Frequency_of_inability_to_cease_drinking_in_last_year",
  "genebass_F1x_Date_F11_first_reported_(mental_and_behavioural_disorders_due_to_use_of_opioids)",
  # F2
  "DisGeNET_F2x_Psychotic_Disorders",
  "DisGeNET_F2x_Schizophrenia",
  "GWASCatalog_F2x_schizophrenia,sex interaction measurement",
  "GWASCatalog_F2x_schizophrenia,intelligence,self reported educational attainment",
  "genebass_F2x_Date_F22_first_reported_(persistent_delusional_disorders)",
  "genebass_F2x_Date_F23_first_reported_(acute_and_transient_psychotic_disorders)",
  # F3
  "DisGeNET_F3x_Bipolar_Disorder",
  "DisGeNET_F3x_Depressive_disorder",
  "GWASCatalog_F3x_age at onset,major depressive disorder",
  "genebass_F3x_Date_F33_first_reported_(recurrent_depressive_disorder)",
  "genebass_F3x_Recent_trouble_concentrating_on_things",
  "genebass_F3x_Recent_feelings_of_inadequacy",
  # F4
  "DisGeNET_F4x_Anxiety_Disorders",
  "DisGeNET_F4x_Panic_Disorder",
  "GWASCatalog_F4x_post-traumatic stress disorder",
  "GWASCatalog_F4x_panic disorder",
  "genebass_F4x_Recent_feelings_of_foreboding",
  "genebass_F4x_Date_F40_first_reported_(phobic_anxiety_disorders)",
  # F5
  "DisGeNET_F5x_Erectile_dysfunction",
  "DisGeNET_F5x_Pediatric_failure_to_thrive",
  "GWASCatalog_F5x_bulimia nervosa",
  "GWASCatalog_F5x_insomnia",
  "genebass_F5x_Date_F53_first_reported_(mental_and_behavioural_disorders_associated_with_the_puerperium,_not_elsewhere_classified)",
  "genebass_F5x_Date_F51_first_reported_(nonorganic_sleep_disorders)",
  # F6
  "GWASCatalog_F6x_childhood gender nonconformity",
  "GWASCatalog_F6x_internet addiction disorder",
  "genebass_F6x_Date_F66_first_reported_(psychological_and_behavioural_disorders_associated_with_sexual_development_and_orientation)",
  "genebass_F6x_Date_F60_first_reported_(specific_personality_disorders)",
  # F7
  "DisGeNET_F7x_Mild_intellectual_disability",
  "DisGeNET_F7x_Intellectual_Disability",
  "genebass_F7x_Fluid_intelligence_score",
  "genebass_F7x_Attempted_fluid_intelligence_(FI)_test.",
  # F8
  "DisGeNET_F8x_Learning_Disorders",     
  "DisGeNET_F8x_Learning_Disabilities",
  "GWASCatalog_F8x_autism spectrum disorder",
  "GWASCatalog_F8x_dyslexia",
  "genebass_F8x_Date_F81_first_reported_(specific_developmental_disorders_of_scholastic_skills)",
  "genebass_F8x_Date_F80_first_reported_(specific_developmental_disorders_of_speech_and_language)",
  # F9
  "DisGeNET_F9x_Stereotypic_Movement_Disorder",
  "DisGeNET_F9x_Gilles_de_la_Tourette_syndrome",
  "GWASCatalog_F9x_attention deficit hyperactivity disorder,risk-taking behaviour",
  "GWASCatalog_F9x_attention deficit hyperactivity disorder"
)

# ##############################################################################
# ---- heatmap ----
# ##############################################################################
GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$list$chi2_value_matrix

GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$list$log_chi2_value_matrix <- log2(GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$list$chi2_value_matrix + 1)

GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$list$log_number_overlap_matrix <- log2(GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$list$number_overlap_matrix + 1)
GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$list$log10_number_overlap_matrix <- log10(GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$list$number_overlap_matrix + 1)


GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$list$sqrt_number_overlap_matrix <- sqrt(GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$list$number_overlap_matrix)

GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$list$log_p_value_matrix <- -log10(GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$list$p_value_matrix) 




GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$list$score_A <-
  GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$list$log_chi2_value_matrix *
  GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$list$log_number_overlap_matrix

GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$list$score_B <-
  GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$list$log_chi2_value_matrix *
  GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$list$sqrt_number_overlap_matrix

GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$list$score_C <-
  GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$list$log_chi2_value_matrix *
  GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$list$number_overlap_matrix

GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$list$score_D <-
  GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$list$log_chi2_value_matrix *
  GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$list$log_number_overlap_matrix *
  GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$list$log_p_value_matrix 

GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$list$score_E <-
  GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$list$log2_odds_ratio_matrix *
  GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$list$log_number_overlap_matrix 

GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$list$score_F <-
  GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$list$log_number_overlap_matrix *
  GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$list$log_p_value_matrix 

GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$list$score_G <-
  GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$list$log10_number_overlap_matrix *
  GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$list$log_p_value_matrix 



patterns <- c("F0x","F1x","F2x","F3x","F4x","F5x","F6x","F7x","F8x","F9x")


GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$df %>% 
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
    icd10_category = map_chr(matched, ~ if(length(.x)==1) .x else NA_character_),
    
    # ---- zamiana nazw ----
    Var2 = recode(Var2,
                  "global_GR_genes_globalUp5TissuesDerivedCells"   = "systemicUp",
                  "global_GR_genes_globalDown5TissuesDerivedCells" = "systemicDown"
    )
  ) %>% 
  select(-c(fdr, fdr_value, matched, n_unique_patterns)) %>% 
  group_by(Var2, icd10_category) %>% 
  slice_max(order_by = combine_score, n = 10, with_ties = FALSE) %>%
  # filter(icd10_category %in% c("F0x", "F1x", "F2x", "F3x", "F4x", "F8x")) %>%
  ungroup() %>% 
  # filter(icd10_category ==  "F2x") %>% 
  # filter(p_value < 0.05) %>% 
  # filter(gene_overlap_count < 3) %>% 
  # filter(icd10_category == "F9x") %>% 
  # .$Var1
  group_by(icd10_category) %>% 
  filter(combine_score > 0) %>% 
  nest() %>%
  mutate(mean_cs = map(data, ~mean(.x$combine_score))) %>%
  unnest(mean_cs)
  select(Var1, icd10_category) %>%
  unique() %>% 
  # head %>% 
  arrange(icd10_category) %>% 
  .$Var1 -> GWASCatalogDisGeNETgenebass_phenotypes
  
  
  
GWASCatalogDisGeNETgenebass_phenotypes <- c(
  # F0x
  "genebass_F0x_Date_F05_first_reported_(delirium,_not_induced_by_alcohol_and_other_psychoactive_substances)",
  "DisGeNET_F0x_Delirium,_Dementia,_Amnestic,_Cognitive_Disorders", 
  "GWASCatalog_F0x_dementia", 
  
  "genebass_F0x_Date_F06_first_reported_(other_mental_disorders_due_to_brain_damage_and_dysfunction_and_to_physical_disease)",
  "genebass_F0x_Date_F01_first_reported_(vascular_dementia)", 
  "genebass_F0x_Date_F00_first_reported_(dementia_in_alzheimer's_disease)",
  
  ###                                           
  "genebass_F0x_Date_F03_first_reported_(unspecified_dementia)",
  "genebass_F0x_Date_F09_first_reported_(unspecified_organic_or_symptomatic_mental_disorder)",                                
  "GWASCatalog_F0x_vascular dementia",
  
  # F1x
  "DisGeNET_F1x_Alcoholic_Intoxication,_Chronic",                                                                             
  "DisGeNET_F1x_Alcoholic_Intoxication", 
  "genebass_F1x_Date_F17_first_reported_(mental_and_behavioural_disorders_due_to_use_of_tobacco)",                            

  "genebass_F1x_Frequency_of_inability_to_cease_drinking_in_last_year",                                                       
  "genebass_F1x_Date_F11_first_reported_(mental_and_behavioural_disorders_due_to_use_of_opioids)",                            
  "genebass_F1x_Age_when_last_took_cannabis",
  ###
  "DisGeNET_F1x_Marijuana_Abuse",                     
  "genebass_F1x_Frequency_of_memory_loss_due_to_drinking_alcohol_in_last_year",  
  
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
 ###                                                                            
 "DisGeNET_F3x_Unipolar_Depression",                                                                                         
 "DisGeNET_F3x_Major_depression,_single_episode_(disorder)",
 
 # F4x
 "genebass_F4x_Date_F48_first_reported_(other_neurotic_disorders)",                                                         
 "genebass_F4x_Date_F41_first_reported_(other_anxiety_disorders)",                                                           
 "genebass_F4x_Date_F44_first_reported_(dissociative_[conversion]_disorders)",
 
 "genebass_F4x_Recent_feelings_of_foreboding",                                                                               
 "genebass_F4x_Date_F40_first_reported_(phobic_anxiety_disorders)" ,                                                         
 "genebass_F4x_Recent_easy_annoyance_or_irritability" ,         
 ###
 "GWASCatalog_F4x_post-traumatic stress disorder" ,                                                                          
 "genebass_F4x_Recent_feelings_or_nervousness_or_anxiety" ,    
 
 # F8x
 "DisGeNET_F8x_Learning_Disorders",                                                                                          
 "DisGeNET_F8x_Learning_Disabilities",                                                                                      
 "DisGeNET_F8x_Autism_Spectrum_Disorders",
 
 "genebass_F8x_Date_F81_first_reported_(specific_developmental_disorders_of_scholastic_skills)",                             
 "DisGeNET_F8x_Developmental_Disabilities",                                                                                  
 "genebass_F8x_Date_F80_first_reported_(specific_developmental_disorders_of_speech_and_language)",                           
 
 ####
 "DisGeNET_F8x_Developmental_delay",                                                                                         
 "DisGeNET_F8x_Autistic_Disorder",                                                                                           
 "DisGeNET_F8x_Neurodevelopmental_Disorders" ,
 
 #
 "genebass_F6x_Date_F66_first_reported_(psychological_and_behavioural_disorders_associated_with_sexual_development_and_orientation)",
 "genebass_F6x_Date_F60_first_reported_(specific_personality_disorders)",
 # F7
 "DisGeNET_F7x_Mild_intellectual_disability",
 "DisGeNET_F7x_Intellectual_Disability",
 "genebass_F7x_Fluid_intelligence_score",
 "genebass_F7x_Attempted_fluid_intelligence_(FI)_test."
)


  
 
heatmap_overlap_log2OR_complex(
  data_list = GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed,
  data_type = "original_data",
  
  # 🔺 triangle mode pełny
  triangle_mode = "full",
  color_key = "score_G",
  # label_key = "same_as_color",
  
  # 🔹 filtrowanie
  rows_to_filter = GWASCatalogDisGeNETgenebass_phenotypes,
  
  row_order_original = GWASCatalogDisGeNETgenebass_phenotypes,
  
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
    "global_GR_genes_globalDown5TissuesDerivedCells" =  "systemicDown",
    "global_GR_genes_globalUp5TissuesDerivedCells" =  "systemiclUp"
  ),
  row_dend_height = unit(20, "mm"),
  col_dend_height = unit(20, "mm"),
  tile_gap = 1,
  # save_to_svg = "results_v2/overlap/overlap_GWASCatalogDisGeNETgenebass/figures/heatmap_GRsystemic_GWASCatalogDisGeNETgenebass_log10pvalue_log10overlap_v2_13.02.2026.svg",
  svg_width = 6.625,
  svg_height = 22,
  row_names_width = unit(90, "mm"),
  col_names_height = unit(50, "mm"),
  force_create_directory = TRUE
)


# ##############################################################################
# ---- save resuolts to XLSX file ----
# ##############################################################################
GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$df %>% 
  select(c(Var1, Var2, p_value, chi2, odds_ratio, log2_odds_ratio, gene_overlap_count, overlap_genes)) %>% 
  
  rename(phenotype = "Var1") %>% 
  rename(grSignature = "Var2") %>% 
  mutate(grSignature = str_replace_all(grSignature, "global_GR_genes_globalUp5TissuesDerivedCells", "systemicUp")) %>% 
  mutate(grSignature = str_replace_all(grSignature, "global_GR_genes_globalDown5TissuesDerivedCells", "systemicDown")) %>% 
  mutate(
    phenotypeSource = case_when(
      # 1) warunek prosty
      grepl("DisGeNET_", phenotype)  ~ "DisGeNET",
      grepl("GWASCatalog_", phenotype)  ~ "GWASCatalog",
      grepl("genebass_", phenotype)  ~ "genebass",
      TRUE ~ "other"
    )
  ) %>% 
  select(phenotypeSource, everything()) %>% 
  mutate(phenotype = str_remove(phenotype, "DisGeNET_")) %>% 
  mutate(phenotype = str_remove(phenotype, "GWASCatalog_")) %>% 
  mutate(phenotype = str_remove(phenotype, "genebass_")) %>% 
  arrange(phenotypeSource)



GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$df %>% 
  select(c(Var1, Var2, p_value, chi2, odds_ratio, log2_odds_ratio, gene_overlap_count, overlap_genes)) %>% 
  
  rename(phenotype = "Var1") %>% 
  rename(grSignature = "Var2") %>% 
  mutate(grSignature = str_replace_all(grSignature, "global_GR_genes_globalUp5TissuesDerivedCells", "systemicUp")) %>% 
  mutate(grSignature = str_replace_all(grSignature, "global_GR_genes_globalDown5TissuesDerivedCells", "systemicDown")) %>% 
  mutate(
    phenotypeSource = case_when(
      # 1) warunek prosty
      grepl("DisGeNET_", phenotype)  ~ "DisGeNET",
      grepl("GWASCatalog_", phenotype)  ~ "GWASCatalog",
      grepl("genebass_", phenotype)  ~ "genebass",
      TRUE ~ "other"
    )
  ) %>% 
  select(phenotypeSource, everything()) %>% 
  mutate(phenotype = str_remove(phenotype, "DisGeNET_")) %>% 
  mutate(phenotype = str_remove(phenotype, "GWASCatalog_")) %>% 
  mutate(phenotype = str_remove(phenotype, "genebass_")) %>% 
  arrange(phenotypeSource) %>% 
  filter(p_value < 0.05 & log2_odds_ratio > 0, gene_overlap_count >= 3)


GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$df %>% 
  select(c(Var1, Var2, p_value, chi2, odds_ratio, log2_odds_ratio, gene_overlap_count, overlap_genes)) %>% 
  
  rename(phenotype = "Var1") %>% 
  rename(grSignature = "Var2") %>% 
  mutate(grSignature = str_replace_all(grSignature, "global_GR_genes_globalUp5TissuesDerivedCells", "systemicUp")) %>% 
  mutate(grSignature = str_replace_all(grSignature, "global_GR_genes_globalDown5TissuesDerivedCells", "systemicDown")) %>% 
  filter(phenotype %in% GWASCatalogDisGeNETgenebass_phenotypes) %>% 
  mutate(
    phenotypeSource = case_when(
      # 1) warunek prosty
      grepl("DisGeNET_", phenotype)  ~ "DisGeNET",
      grepl("GWASCatalog_", phenotype)  ~ "GWASCatalog",
      grepl("genebass_", phenotype)  ~ "genebass",
      TRUE ~ "other"
    )
  ) %>% 
  select(phenotypeSource, everything()) %>% 
  mutate(phenotype = str_remove(phenotype, "DisGeNET_")) %>% 
  mutate(phenotype = str_remove(phenotype, "GWASCatalog_")) %>% 
  mutate(phenotype = str_remove(phenotype, "genebass_")) %>% 
  arrange(phenotypeSource)




# ##############################################################################
suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(openxlsx)
})

out_dir <- "/home/mateusz/projects/ifpan-michkor-gr/results_v2/overlap/overlap_GWASCatalogDisGeNETgenebass/tables"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_file <- file.path(out_dir, "file_GRsystemic_GWASCatalogDisGeNETgenebass_log2OR_30.01.2026.xlsx")

# =========================
# 1) allResults
# =========================
allResults <- GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$df %>%
  select(c(Var1, Var2, p_value, chi2, odds_ratio, log2_odds_ratio, gene_overlap_count, overlap_genes)) %>%
  rename(phenotype = "Var1") %>%
  rename(grSignature = "Var2") %>%
  mutate(grSignature = str_replace_all(grSignature, "global_GR_genes_globalUp5TissuesDerivedCells", "systemicUp")) %>%
  mutate(grSignature = str_replace_all(grSignature, "global_GR_genes_globalDown5TissuesDerivedCells", "systemicDown")) %>%
  mutate(
    phenotypeSource = case_when(
      grepl("DisGeNET_", phenotype) ~ "DisGeNET",
      grepl("GWASCatalog_", phenotype) ~ "GWASCatalog",
      grepl("genebass_", phenotype) ~ "genebass",
      TRUE ~ "other"
    )
  ) %>%
  select(phenotypeSource, everything()) %>%
  mutate(phenotype = str_remove(phenotype, "DisGeNET_")) %>%
  mutate(phenotype = str_remove(phenotype, "GWASCatalog_")) %>%
  mutate(phenotype = str_remove(phenotype, "genebass_")) %>%
  arrange(phenotypeSource)

# =========================
# 2) significantResults
# =========================
significantResults <- GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$df %>%
  select(c(Var1, Var2, p_value, chi2, odds_ratio, log2_odds_ratio, gene_overlap_count, overlap_genes)) %>%
  rename(phenotype = "Var1") %>%
  rename(grSignature = "Var2") %>%
  mutate(grSignature = str_replace_all(grSignature, "global_GR_genes_globalUp5TissuesDerivedCells", "systemicUp")) %>%
  mutate(grSignature = str_replace_all(grSignature, "global_GR_genes_globalDown5TissuesDerivedCells", "systemicDown")) %>%
  mutate(
    phenotypeSource = case_when(
      grepl("DisGeNET_", phenotype) ~ "DisGeNET",
      grepl("GWASCatalog_", phenotype) ~ "GWASCatalog",
      grepl("genebass_", phenotype) ~ "genebass",
      TRUE ~ "other"
    )
  ) %>%
  select(phenotypeSource, everything()) %>%
  mutate(phenotype = str_remove(phenotype, "DisGeNET_")) %>%
  mutate(phenotype = str_remove(phenotype, "GWASCatalog_")) %>%
  mutate(phenotype = str_remove(phenotype, "genebass_")) %>%
  arrange(phenotypeSource) %>%
  filter(p_value < 0.05 & log2_odds_ratio > 0, gene_overlap_count >= 3)

# =========================
# 3) heatmapResults
# =========================
heatmapResults <- GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$df %>%
  select(c(Var1, Var2, p_value, chi2, odds_ratio, log2_odds_ratio, gene_overlap_count, overlap_genes)) %>%
  rename(phenotype = "Var1") %>%
  rename(grSignature = "Var2") %>%
  mutate(grSignature = str_replace_all(grSignature, "global_GR_genes_globalUp5TissuesDerivedCells", "systemicUp")) %>%
  mutate(grSignature = str_replace_all(grSignature, "global_GR_genes_globalDown5TissuesDerivedCells", "systemicDown")) %>%
  filter(phenotype %in% GWASCatalogDisGeNETgenebass_phenotypes) %>%
  mutate(
    phenotypeSource = case_when(
      grepl("DisGeNET_", phenotype) ~ "DisGeNET",
      grepl("GWASCatalog_", phenotype) ~ "GWASCatalog",
      grepl("genebass_", phenotype) ~ "genebass",
      TRUE ~ "other"
    )
  ) %>%
  select(phenotypeSource, everything()) %>%
  mutate(phenotype = str_remove(phenotype, "DisGeNET_")) %>%
  mutate(phenotype = str_remove(phenotype, "GWASCatalog_")) %>%
  mutate(phenotype = str_remove(phenotype, "genebass_")) %>%
  arrange(phenotypeSource)

# =========================
# Zapis do jednego XLSX (3 sheety) + freeze top row
# =========================
wb <- createWorkbook()

addWorksheet(wb, "allResults")
writeData(wb, "allResults", allResults)
freezePane(wb, "allResults", firstRow = TRUE)

addWorksheet(wb, "significantResults")
writeData(wb, "significantResults", significantResults)
freezePane(wb, "significantResults", firstRow = TRUE)

addWorksheet(wb, "heatmapResults")
writeData(wb, "heatmapResults", heatmapResults)
freezePane(wb, "heatmapResults", firstRow = TRUE)

saveWorkbook(wb, out_file, overwrite = TRUE)

message("Saved: ", out_file)



