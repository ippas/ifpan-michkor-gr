
#
# ##############################################################################
# ---- uses data ----
# ##############################################################################
grSignatures_6cellTypes_edgerCPM50p0.01log2FC1 <- readRDS("results_v2/GR_signatures/grSignatures_6cellTypes_edgerCPM50p0.01log2FC1.rds")

disgenet_mentalDisorders <- readRDS("data/databases/disgenet/disgenet_mentalDisordersF03.rds")
# disgenet_mentalDisorders <- readRDS("data/databases/disgenet/disgenet_MSH_cardiovascularDiseasesC14.rds")

random_8geneLists <-  generate_random_gene_lists(hgnc_symbols_vector_v110, n_lists = 8, length = 50)
# ##############################################################################
# ---- functions ----
# ##############################################################################
filter_min_vector_length <- function(x, min_len = 10) {
  Filter(function(v) length(v) >= min_len, x)
}


# ##############################################################################
# ---- chi2 ----
# ##############################################################################
disgenetMentalHealth_GrSignatures_overlapChi2 <- run_full_overlap_analysis(
  gene_lists = c(grSignatures_6cellTypes_edgerCPM50p0.01log2FC1, 
                 random_8geneLists,
                 disgenet_mentalDisorders$geneLists_scoreMin0.5 %>%
                   filter_min_vector_length(min_len = 10)
  ),
  total_genes = hgnc_symbols_vector_v110,
  cols_to_filter = c(names(grSignatures_6cellTypes_edgerCPM50p0.01log2FC1), names(random_8geneLists)),
  rows_to_filter = names(disgenet_mentalDisorders$geneLists_scoreMin0.5 %>%
                           filter_min_vector_length(min_len = 10)),
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
  filter(Var2 %in% names(grSignatures_6cellTypes_edgerCPM50p0.01log2FC1)) %>%
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
  # cols_to_filter = disgenetMentalHealth_grSignatures_vector_p0.05,
  
  # 🎨 skala kolorów
  color_scale_range = c(-3, 3),
  text_contrast_range = c(-30, 2.5),
  palette = c("#c6d3e3", "white", "darkred"),
  
  # 🔲 prostokąty istotności
  color_rects = c("#97C426", "#2F4603"),
  
  # 📊 klastrowanie
  cluster_rows = T,
  cluster_cols = F,
  show_dendrograms = TRUE,
  rect_lwd = 2.5,
  p_thresholds = c(0.05, 0.01),
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
  # save_to_svg = "results_v2/overlap/pgc_overlap/figures/heatmap_SignifGrSignaturesPGC_log2OR.svg",
  # save_to_svg = "results_v2/overlap/disgenet_overlap/figures/heatmap_AllGrSignatures_disgenetMentalHealthMinScore0_log2OR.svg",
  svg_width = 10.5,
  svg_height = 12,
  row_names_width = unit(90, "mm"),
  col_names_height = unit(90, "mm"),
  force_create_directory = TRUE
)


# ##############################################################################
disgenetMentalHealth_GrSignatures_overlapChi2$processed$original_data$df %>% 
  filter(gene_overlap_count > 2) %>% 
  filter(p_value < 0.05) %>% 
  filter(Var2 %in% sig_names) %>%
  group_by(Var2) %>%
  nest() %>% 
  mutate(n_association = map(data, ~ .x %>% nrow)) %>% 
  mutate(all_genes = map(data, ~ .x %>% .$overlap_genes %>% strsplit(",") %>% unlist() %>% unique %>% paste(., collapse = "|"))) %>% 
  mutate(n_genes = map(data, ~ .x %>% .$overlap_genes %>% strsplit(",") %>% unlist() %>% unique %>% length)) %>% 
  unnest(n_association, n_genes) %>% 
  select(c(Var2, data, n_association, n_genes, all_genes)) %>% 
  select(-data) %>% 
  mutate(
    Var2 = case_when(
      str_detect(Var2, "^global_GR_genes_globalUp")   ~ "globalUp",
      str_detect(Var2, "^global_GR_genes_globalDown") ~ "globalDown",
      TRUE ~ str_remove(Var2, "^minusGlobalUpDown5TissuesDerivedCells_")
    )
  ) 


disgenetMentalHealth_GrSignatures_overlapChi2$processed$original_data$df %>% 
  filter(gene_overlap_count > 2) %>% 
  filter(p_value < 0.05) %>% 
  filter(grepl("FKBP5", overlap_genes))



output <- summarise_disgenet_multiscore(
  score_values = disgenet_score_min,
  flat_signatures = flat_allGrSignatures_31.10.2025,
  sig_names = sig_names,
  disgenet_object = disgenet_mentalDisorders,
  total_genes = hgnc_symbols_vector_v110
)

output$summary_long_chi2 %>% 
  filter(score_min == 0.5) %>% t


# raw_list <- disgenet_multiDisease2genes(
#   disease_ids  = c("C0036341"),
#   database     = "CLINVAR",
#   score        = c(0,1),
#   verbose      = TRUE
# )
# raw_list$Schizophrenia %>% dim
