# ##############################################################################
# ---- reada genebass data ----
# ##############################################################################
# genebass_MentalAndBehaviouralDisorders <- read_tsv("data/genebass_v2/all_categories_SKAT/genebass_SKAT_Health-related_outcomes_-_First_occurrences_-_Mental_and_behavioural_disorders.tsv.bgz") %>%
# genebass_MentalAndBehaviouralDisorders <- read_tsv("data/genebass_v2/all_categories_SKATO/genebass_SKATO_Health-related_outcomes_-_First_occurrences_-_Mental_and_behavioural_disorders.tsv.bgz") %>% 
genebass_MentalAndBehaviouralDisorders <- read_tsv("data/genebass_v2/all_categories_burden/genebass_burden_Health-related_outcomes_-_First_occurrences_-_Mental_and_behavioural_disorders.tsv.bgz") %>%
# genebass_MentalAndBehaviouralDisorders <- read_tsv("data/genebass_v2/all_categories_SKAT/genebass_SKAT_UK_Biobank_Assessment_Centre_-_Touchscreen_-_Psychosocial_factors_-_Mental_health.tsv.bgz") %>%
# genebass_MentalAndBehaviouralDisorders <- read_tsv("data/genebass_v2/all_categories_burden/genebass_burden_UK_Biobank_Assessment_Centre_-_Touchscreen_-_Psychosocial_factors_-_Mental_health.tsv.bgz") %>%
# genebass_MentalAndBehaviouralDisorders <- read_tsv("data/genebass_v2/all_categories_burden/genebass_burden_Health-related_outcomes_-_First_occurrences_-_Nervous_system_disorders.tsv.bgz") %>%
  select(-c(pvalue_test, pvalue_threshold, heritability)) %>% 
  filter(annotation == "pLoF") %>%
  filter(pvalue < 0.05) %>% 
  mutate(description_format = str_replace_all(description, " ", "_"))

genebass_MentalAndBehaviouralDisorders$description %>% unique %>% length()

gene_list <- genebass_MentalAndBehaviouralDisorders %>% 
  filter(pvalue < 0.05) %>% 
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
  group_by(description_format) %>%
  summarise(genes = list(unique(gene_symbol)), .groups = "drop") %>%
  deframe()

random_8geneLists <-  generate_random_gene_lists(hgnc_symbols_vector_v110, n_lists = 8, length = 150)
# ============================================================
# 🧩 2️⃣ Uruchomienie analizy overlap (dla nowych sygnatur)
# ============================================================

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

genebasGrSignatures_overlapChi2 <- run_full_overlap_analysis(
  gene_lists = c(flat_allGrSignatures_18.11.2025[sig_names], gene_list, random_8geneLists),
  total_genes = hgnc_symbols_vector_v110,
  cols_to_filter = c(sig_names, names(random_8geneLists)),
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
  filter(p_value < 0.05) %>% 
  filter(gene_overlap_count  >= 3) %>%
  filter(log2_odds_ratio > 0) %>% 
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
  # cols_to_filter = pgc_grSignatures_vector_p0.05,
  
  # 🎨 skala kolorów
  color_scale_range = c(-3, 3),
  text_contrast_range = c(-30, 4.9),
  palette = c("#c6d3e3", "white", "darkred"),
  
  # 🔲 prostokąty istotności
  color_rects = c("#97C426", "#2F4603"),
  
  # 📊 klastrowanie
  cluster_rows = T,
  cluster_cols = F,
  show_dendrograms = TRUE,
  rect_lwd = 2.5,
  col_mapper = c(
    "minusGlobalUpDown5TissuesDerivedCells_BloodCellsUp" = "BloodCellsUp",
    "minusGlobalUpDown5TissuesDerivedCells_LungCellsUp" = "LungCellsUp",
    "minusGlobalUpDown5TissuesDerivedCells_NeuralCellsUp" = "NeuralCellsUp",
    "minusGlobalUpDown5TissuesDerivedCells_BloodCellsDown" = "BloodCellsDown",
    "minusGlobalUpDown5TissuesDerivedCells_LungCellsDown" = "LungCellsDown",
    "minusGlobalUpDown5TissuesDerivedCells_NeuralCellsDown" = "NeuralCellsDown",
    "global_GR_genes_globalDown5TissuesDerivedCells" =  "globalDown",
    "global_GR_genes_globalUp5TissuesDerivedCells" =  "globalUp"
  ),
  row_dend_height = unit(20, "mm"),
  col_dend_height = unit(20, "mm"),
  tile_gap = 1,
  # save_to_svg = "results_v2/overlap/pgc_overlap/figures/heatmap_allGrSignaturesPGC_log2OR.svg",
  svg_width = 10.5, 
  svg_height = 10,
  row_names_width = unit(90, "mm"),
  col_names_height = unit(50, "mm"),
  force_create_directory = TRUE
)



  
