input_dir <- "data/genebass_v2/all_categories_SKAT/"

files <- c(
  list.files(
    path = input_dir,
    pattern = "SKAT_Online_follow-up_.*Mental_health.*\\.tsv\\.bgz$",
    full.names = TRUE
  ),
  file.path(input_dir, "genebass_SKAT_Health-related_outcomes_-_First_occurrences_-_Mental_and_behavioural_disorders.tsv.bgz"),
  file.path(input_dir, "genebass_SKAT_Health-related_outcomes_-_First_occurrences_-_Nervous_system_disorders.tsv.bgz")
)

files <- c(
  list.files(
    path = input_dir,
    pattern = "SKAT_Online_follow-up_.*Mental_health.*\\.tsv\\.bgz$",
    full.names = TRUE
  )
)

mapperPhenotyeps2icd10 <- read_excel_sheets("data/mapperPhenotypesICD10/manualMapping_phenotype2ICD10categories_24.01.2026.xlsx") %>% 
  set_names(c( "GWASCatalog", "DisGeNET","genebass"))

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

mapperPhenotyeps2icd10$genebass %>% 
  mutate(phenotype = str_replace_all(phenotype, " ", "_")) %>% 
  mutate(names = paste0(subcategory, "_", phenotype)) %>% 
  mutate(names = str_remove(names, "NA_")) %>% 
  mutate(names = paste0(mapping2ICD10, "_", names)) %>% 
  select(c(phenotype, names)) -> mapperPhenotypesGenebass


# ##############################################################################
# ---- analysis ----
# ##############################################################################

gene_list <- genebass_online_mental %>% 
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
  group_by(names) %>%
  summarise(genes = list(unique(gene_symbol)), .groups = "drop") %>%
  deframe()

random_8geneLists <-  generate_random_gene_lists(hgnc_symbols_vector_v110, n_lists = 8, length = 150)
# ============================================================
# 🧩 2️⃣ Uruchomienie analizy overlap (dla nowych sygnatur)
# ============================================================

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



genebasGrSignatures_overlapChi2 <- run_full_overlap_analysis(
  # gene_lists = c(flat_allGrSignatures_18.11.2025[sig_names], gene_list, random_8geneLists),
  gene_lists = c(flat_allGrSignatures_18.11.2025[sig_names], gene_list),
  total_genes = hgnc_symbols_vector_v110,
  cols_to_filter = c(sig_names),
  rows_to_filter = names(gene_list),
  plot_title_or = "",
  triangle_mode = "full",
  fdr_threshold = 1,
  data_type = "original_data",
  verbose = FALSE
)


# ============================================================
# 🧠 4️⃣ Wyznacz listę nazw sygnatur z istotnymi overlapami
# ============================================================

genebasGrSignatures_overlapChi2$processed$original_data$df %>% 
  filter(p_value < 0.05) %>% 
  filter(gene_overlap_count  >= 3) %>% 
  filter(log2_odds_ratio > 0) %>% 
  .$Var2 %>% as.character() %>% unique() -> genebass_grSignatures_vector_p0.05

genebasGrSignatures_overlapChi2$processed$original_data$df %>% 
  filter(grepl("F1x|F2x|F3x|F4x", Var1)) %>% 
  # filter(grepl("F4x", Var1)) %>% 
  filter((p_value < 0.05 & gene_overlap_count >= 3 & log2_odds_ratio > 0) | grepl("F2x", Var1)) %>%
  # filter(gene_overlap_count >= 3) %>%
  # filter(log2_odds_ratio > 0) %>%
  .$Var1 %>% as.character() %>% unique()  -> genebass_phenotypes_vector_p0.05

# ##############################################################################
# ---- complex heatmap ----
# ##############################################################################

# ---- All GrSignatures ----
heatmap_overlap_log2OR_complex(
  data_list = genebasGrSignatures_overlapChi2$processed,
  data_type = "original_data",
  
  # 🔺 triangle mode pełny
  triangle_mode = "full",
  
  # 🔹 filtrowanie
  rows_to_filter = genebass_phenotypes_vector_p0.05,
  
  col_order_original = c("globalUp", "globalDown" 
                         # "NeuralCellsUp", "NeuralCellsDown",
                         # "BloodCellsUp", "BloodCellsDown", "LungCellsUp", "LungCellsDown"
                         ),
  
  row_order_original = c(
    "F1x_Alcohol use_Frequency_of_inability_to_cease_drinking_in_last_year",
    "F1x_Cannabis use_Age_when_last_took_cannabis",
    "F3x_Depression_Recent_feelings_of_inadequacy",
    "F3x_Depression_Recent_trouble_concentrating_on_things",
    "F4x_Anxiety_Recent_feelings_of_foreboding",
    "F2x_Unusual and psychotic experiences_Ever_heard_an_un-real_voice",
    "F2x_Unusual and psychotic experiences_Ever_believed_in_un-real_communications_or_signs",
    "F2x_Unusual and psychotic experiences_Ever_believed_in_an_un-real_conspiracy_against_self",
    "F2x_Unusual and psychotic experiences_Ever_seen_an_un-real_vision",
    "F2x_Unusual and psychotic experiences_Ever_prescribed_a_medication_for_unusual_or_psychotic_experiences",
    "F2x_Unusual and psychotic experiences_Ever_talked_to_a_health_professional_about_unusual_or_psychotic_experiences"
  ),
  # cols_to_filter = pgc_grSignatures_vector_p0.05,
  
  # 🎨 skala kolorów
  color_scale_range = c(-5, 5),
  text_contrast_range = c(-30, 4.9),
  palette = c("#c6d3e3", "white", "darkred"),
  
  # 🔲 prostokąty istotności
  color_rects = c("#97C426", "#2F4603"),
  
  # 📊 klastrowanie
  cluster_rows = F,
  cluster_cols = F,
  show_dendrograms = TRUE,
  rect_lwd = 2.5,
  col_mapper = c(
    # "minusGlobalUpDown5TissuesDerivedCells_BloodCellsUp" = "BloodCellsUp",
    # "minusGlobalUpDown5TissuesDerivedCells_LungCellsUp" = "LungCellsUp",
    # "minusGlobalUpDown5TissuesDerivedCells_NeuralCellsUp" = "NeuralCellsUp",
    # "minusGlobalUpDown5TissuesDerivedCells_BloodCellsDown" = "BloodCellsDown",
    # "minusGlobalUpDown5TissuesDerivedCells_LungCellsDown" = "LungCellsDown",
    # "minusGlobalUpDown5TissuesDerivedCells_NeuralCellsDown" = "NeuralCellsDown",
    "global_GR_genes_globalDown5TissuesDerivedCells" =  "globalDown",
    "global_GR_genes_globalUp5TissuesDerivedCells" =  "globalUp"
  ),
  row_dend_height = unit(20, "mm"),
  col_dend_height = unit(20, "mm"),
  tile_gap = 1,
  save_to_svg = "results_v2/overlap/genebass_overlap/figures/heatmap_systemicGrSignaturesGenebass_log2OR_25.01.2026.svg",
  svg_width = 6.62, 
  svg_height = 7.25,
  row_names_width = unit(90, "mm"),
  col_names_height = unit(50, "mm"),
  force_create_directory = TRUE
)


genebasGrSignatures_overlapChi2$processed$original_data$df %>% 
  filter(p_value < 0.05) %>% 
  filter(gene_overlap_count  >= 3) %>%
  filter(log2_odds_ratio > 0) %>% 
  select(-c(sig_1, fdr, ._key)) %>% 
  mutate(Var2 = str_replace_all(Var2, "minusGlobalUpDown5TissuesDerivedCells_", "")) %>% 
  mutate(Var2 = str_replace_all(Var2, "global_GR_genes_globalDown5TissuesDerivedCells", "globalDown")) %>% 
  write.xlsx("results_v2/overlap/genebass_overlap/tables/genebasMentalHealthFollowUp_GrSignatures_overlapChi2_13.01.2026.xlsx", 
             rowNames = FALSE)



# =========================
# Sheet 1: FILTERED (main)
# =========================
df_filtered <- genebasGrSignatures_overlapChi2$processed$original_data$df %>% 
  filter(p_value < 0.05) %>% 
  filter(gene_overlap_count >= 3) %>%
  filter(log2_odds_ratio > 0) %>% 
  select(-c(fdr, fdr_value)) %>% 
  mutate(
    Var2 = str_replace_all(Var2, "minusGlobalUpDown5TissuesDerivedCells_", ""),
    Var2 = str_replace_all(Var2, "global_GR_genes_globalDown5TissuesDerivedCells", "systemicDown")
  )

# =========================
# Sheet 2: UNFILTERED (reference)
# =========================
df_unfiltered <- genebasGrSignatures_overlapChi2$processed$original_data$df %>% 
  select(-c( fdr, fdr_value )) %>% 
  mutate(
    Var2 = str_replace_all(Var2, "minusGlobalUpDown5TissuesDerivedCells_", ""),
    Var2 = str_replace_all(Var2, "global_GR_genes_globalDown5TissuesDerivedCells", "systemicDown")
  )

# =========================
# Write XLSX
# =========================
wb <- createWorkbook()

addWorksheet(wb, "Filtered_results")
writeData(wb, "Filtered_results", df_filtered)
freezePane(wb, "Filtered_results", firstRow = TRUE)

addWorksheet(wb, "Unfiltered_reference")
writeData(wb, "Unfiltered_reference", df_unfiltered)
freezePane(wb, "Unfiltered_reference", firstRow = TRUE)

saveWorkbook(
  wb,
  file = "results_v2/overlap/genebass_overlap/tables/genebasMentalHealthFollowUp_GrSystemic_overlapChi2_25.01.2026.xlsx",
  overwrite = TRUE
)


genebass_online_mental %>% 
  select(description, category) %>%
  mutate(description = str_replace_all(description, " ", "_")) %>% 
  filter(description %in% genebass_phenotypes_vector_p0.05) %>% 
  unique

# ##############################################################################
# ---- save to RData ----
# ##############################################################################
save(
  hgnc_symbols_vector_v110,
  flat_allGrSignatures_31.10.2025,
  sig_names,
  gene_list,
  file = "data/run_RData_inputs/permutation_fdr/rdata2permuationFDR_genebassMentalHealth.RData"
)
