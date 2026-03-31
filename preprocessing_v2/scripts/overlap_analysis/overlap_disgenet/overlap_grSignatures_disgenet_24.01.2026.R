
#
# ##############################################################################
# ---- uses data ----
# ##############################################################################
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

mapperPhenotyeps2icd10 <- read_excel_sheets("data/mapperPhenotypesICD10/manualMapping_phenotype2ICD10categories_24.01.2026.xlsx") %>% 
  set_names(c( "GWASCatalog", "DisGeNET","genebass"))

mapperPhenotyeps2icd10$DisGeNET %>% 
  mutate(phenotype = str_replace_all(phenotype, " ", "_")) %>%
  mutate(
    mapping2ICD10 = mapping2ICD10 %>%
      str_replace_all("[ \\+\\(\\);]", "_") %>%  # spacja, +, (, ), ;
      str_replace_all("_+", "_") %>%             # wiele _ → jedno _
      str_replace_all("^_|_$", "")               # usuń _ na początku/końcu
  ) %>% 
  mutate(names = paste0(mapping2ICD10, "_", phenotype)) %>% 
  select(c(phenotype, names)) -> mapperPhenotypesDisGeNET


flat_allGrSignatures_31.10.2025

disgenet_mentalDisorders <- readRDS("data/databases/disgenet/disgenet_mentalDisordersF03.rds")
# disgenet_mentalDisorders <- readRDS("data/databases/disgenet/disgenet_MSH_cardiovascularDiseasesC14.rds")

random_8geneLists <-  generate_random_gene_lists(hgnc_symbols_vector_v110, n_lists = 8, length = 150)
# ##############################################################################
# ---- functions ----
# ##############################################################################
filter_min_vector_length <- function(x, min_len = 10) {
  Filter(function(v) length(v) >= min_len, x)
}


gene_list <- disgenet_mentalDisorders$geneLists_scoreMin0 %>%
  filter_min_vector_length(min_len = 10)

gene_list %>% stack %>% 
  set_colnames(c("gene_symbol", "phenotype")) %>% 
  left_join(., mapperPhenotypesDisGeNET, by = "phenotype") %>% 
  select(names, gene_symbol) %>% 
  group_by(names) %>%
  summarise(gene_symbol = list(gene_symbol), .groups = "drop") %>%
  deframe() -> gene_list

# ##############################################################################
# ---- chi2 ----
# ##############################################################################
disgenetMentalHealth_GrSignatures_overlapChi2 <- run_full_overlap_analysis(
  gene_lists = c(flat_allGrSignatures_31.10.2025[sig_names], 
                 random_8geneLists,
                 gene_list
  ),
  total_genes = hgnc_symbols_vector_v110,
  cols_to_filter = c(sig_names, names(random_8geneLists)),
  rows_to_filter = names(gene_list),
  plot_title_or = "",
  triangle_mode = "full",
  fdr_threshold = 1,
  data_type = "original_data",
  verbose = FALSE,
  palette_or = c("#c6d3e3", "white", "darkred"),
  text_contrast_range_or = c(-30, 4.9),
)


# ============================================================
# 📈 5️⃣ Liczba unikalnych genów GR-zależnych w istotnych overlapach
# ============================================================

disgenetMentalHealth_GrSignatures_overlapChi2$processed$original_data$df %>%
  filter(Var2 %in% sig_names) %>%
  group_by(Var2) %>%
  nest() %>%
  mutate(
    data = map(data, ~ .x %>%
                 mutate(fdr = p.adjust(p_value, method = "fdr")))
  ) %>%
  unnest(data) %>%
  filter(gene_overlap_count > 2) %>%
  filter(log2_odds_ratio > 0) %>%
  filter(p_value < 0.05) -> disgenetMentalHealth_GrSignatures_filtered_df

# ============================================================
# 🧠 4️⃣ Wyznacz listę nazw sygnatur z istotnymi overlapami
# ============================================================

disgenetMentalHealth_GrSignatures_filtered_df$Var2 %>% unique() -> disgenetMentalHealth_grSignatures_vector_p0.05

disgenetMentalHealth_GrSignatures_filtered_df$Var1 %>% unique() -> disgenetMentalHealth_phenotypes_vector_p0.05
disgenetMentalHealth_GrSignatures_overlapChi2$processed$original_data$df$Var1 %>% unique -> disgenetMentalHealth_phenotypes_vector_p0.05

disgenetMentalHealth_GrSignatures_overlapChi2$processed$original_data$df %>% 
  filter((p_value < 0.05 & gene_overlap_count >= 3 & log2_odds_ratio > 0) | grepl("F7x", Var1)) %>% 
  filter(!grepl("randomGeneList", Var2)) %>% 
  .$Var1 %>% 
  unique -> disgenetMentalHealth_phenotypes_vector_p0.05
# ##############################################################################
# ---- heatmap ----
# ##############################################################################

heatmap_overlap_log2OR_complex(
  data_list = disgenetMentalHealth_GrSignatures_overlapChi2$processed,
  data_type = "original_data",
  
  # 🔺 triangle mode pełny
  triangle_mode = "full",
  
  # 🔹 filtrowanie
  rows_to_filter = disgenetMentalHealth_phenotypes_vector_p0.05,
  cols_to_filter = disgenetMentalHealth_grSignatures_vector_p0.05,
  
  # row_order_original = c(
  #   "F1x_Alcohol_abuse",
  #   "F1x_Alcoholic_Intoxication,_Chronic",
  #   "F2x_Nonorganic_psychosis",
  #   "F2x_Psychotic_Disorders",
  #   "F2x_Schizophrenia",
  #   "F3x_Bipolar_Disorder",
  #   "F3x_Depressive_disorder",
  #   "F3x_Major_depression,_single_episode_(disorder)",
  #   "F3x_Major_Depressive_Disorder",
  #   "F3x_Mood_Disorders",
  #   "F3x_Unipolar_Depression",
  #   "F8x_Autistic_Disorder",
  #   "G3x_F0x_Alzheimer's_Disease",
  #   "G3x_Huntington_Disease",
  #   "F7x_Intellectual_Disability",
  #   "F7x_Mental_Retardation",
  #   "F7x_Non-specific_syndromic_intellectual_disability",
  #   "F7x_Severe_intellectual_disability_(disorder)"
  # ),
  
  # 🎨 skala kolorów
  color_scale_range = c(-5, 5),
  text_contrast_range = c(-30, 4),
  palette = c("#c6d3e3", "white", "darkred"),
  
  # 🔲 prostokąty istotności
  color_rects = c("#97C426", "#2F4603"),
  
  # 📊 klastrowanie
  cluster_rows = F,
  cluster_cols = F,
  show_dendrograms = TRUE,
  rect_lwd = 2.5,
  p_thresholds = c(0.05, 0.01),
  col_mapper = c(
    # "minusGlobalUpDown5TissuesDerivedCells_BloodCellsUp" = "BloodCellsUp",
    # "minusGlobalUpDown5TissuesDerivedCells_LungCellsUp" = "LungCellsUp",
    # "minusGlobalUpDown5TissuesDerivedCells_NeuralCellsUp" = "NeuralCellsUp",
    # "minusGlobalUpDown5TissuesDerivedCells_BloodCellsDown" = "BloodCellsDown",
    # "minusGlobalUpDown5TissuesDerivedCells_LungCellsDown" = "LungCellsDown",
    # "minusGlobalUpDown5TissuesDerivedCells_NeuralCellsDown" = "NeuralCellsDown",
    "global_GR_genes_globalDown5TissuesDerivedCells" =  "systemicDown",
    "global_GR_genes_globalUp5TissuesDerivedCells" =  "systemiclUp"
  ),
  row_dend_height = unit(20, "mm"),
  col_dend_height = unit(20, "mm"),
  tile_gap = 1,
  # save_to_svg = "results_v2/overlap/pgc_overlap/figures/heatmap_SignifGrSignaturesPGC_log2OR.svg",
  # save_to_svg = "results_v2/overlap/disgenet_overlap/figures/heatmap_systemicGrSignatures_disgenetMentalHealthMinScore0.5_log2OR_25.01.2026.svg",
  svg_width = 6.53,
  svg_height = 9.85,
  row_names_width = unit(90, "mm"),
  col_names_height = unit(50, "mm"),
  force_create_directory = TRUE
)


# =========================
# Sheet 1: FILTERED (main)
# =========================
df_filtered <- disgenetMentalHealth_GrSignatures_overlapChi2$processed$original_data$df %>%
  filter(p_value < 0.05) %>%
  filter(gene_overlap_count >= 3) %>%
  filter(log2_odds_ratio > 0) %>%
  select(-c(fdr, fdr_value)) %>%
  mutate(
    Var2 = str_replace_all(Var2, "minusGlobalUpDown5TissuesDerivedCells_", ""),
    Var2 = str_replace_all(Var2, "global_GR_genes_globalDown5TissuesDerivedCells", "systemicDown"),
    Var2 = str_replace_all(Var2, "global_GR_genes_globalUp5TissuesDerivedCells", "systemicUp")
  )

# =========================
# Sheet 2: UNFILTERED (reference)
# =========================
df_unfiltered <- disgenetMentalHealth_GrSignatures_overlapChi2$processed$original_data$df %>%
  select(-c(fdr, fdr_value)) %>%
  mutate(
    Var2 = str_replace_all(Var2, "minusGlobalUpDown5TissuesDerivedCells_", ""),
    Var2 = str_replace_all(Var2, "global_GR_genes_globalDown5TissuesDerivedCells", "systemicDown"),
    Var2 = str_replace_all(Var2, "global_GR_genes_globalUp5TissuesDerivedCells", "systemicUp")
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
  file = "results_v2/overlap/disgenet_overlap/tables/disgenetMentalHealth_GrSystemic_overlapChi2_25.01.2026.xlsx",
  overwrite = TRUE
)






