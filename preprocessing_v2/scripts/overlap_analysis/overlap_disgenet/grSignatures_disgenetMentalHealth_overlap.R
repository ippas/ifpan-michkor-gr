# ##############################################################################
# ---- uses data ----
# ##############################################################################
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

flat_allGrSignatures_31.10.2025

disgenet_mentalDisorders <- readRDS("data/databases/disgenet/disgenet_mentalDisordersF03.rds")
# disgenet_mentalDisorders <- readRDS("data/databases/disgenet/disgenet_MSH_cardiovascularDiseasesC14.rds")


# ##############################################################################
# ---- functions ----
# ##############################################################################
filter_min_vector_length <- function(x, min_len = 10) {
  Filter(function(v) length(v) >= min_len, x)
}
summarise_disgenet_multiscore <- function(
    score_values,
    flat_signatures,
    sig_names,
    disgenet_object,
    total_genes,
    min_len = 10,
    verbose = TRUE
) {
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tidyr)
  
  signature_order <- c(
    "globalUp", "globalDown",
    "NeuralCellsUp", "NeuralCellsDown",
    "BloodCellsUp", "BloodCellsDown",
    "LungCellsUp", "LungCellsDown"
  )
  
  get_unique_genes <- function(x) {
    x %>%
      strsplit(",") %>%
      unlist() %>%
      trimws() %>%
      unique()
  }
  
  summarise_overlap <- function(df) {
    df %>%
      mutate(
        Var2 = case_when(
          str_detect(Var2, "^global_GR_genes_globalUp")   ~ "globalUp",
          str_detect(Var2, "^global_GR_genes_globalDown") ~ "globalDown",
          TRUE ~ str_remove(Var2, "^minusGlobalUpDown5TissuesDerivedCells_")
        )
      ) %>%
      group_by(Var2) %>%
      summarise(
        n_assoc_p005 = sum(p_value < 0.05 & gene_overlap_count > 2),
        n_assoc_p001 = sum(p_value < 0.01 & gene_overlap_count > 2),
        
        n_genes_p005 = get_unique_genes(overlap_genes[p_value < 0.05 & gene_overlap_count > 2]) %>% length(),
        n_genes_p001 = get_unique_genes(overlap_genes[p_value < 0.01 & gene_overlap_count > 2]) %>% length(),
        n_genes_all  = get_unique_genes(overlap_genes) %>% length(),
        
        .groups = "drop"
      )
  }
  
  raw_chi2_results <- list()
  summary_by_score_chi2 <- list()
  
  for (sc in score_values) {
    
    if (verbose) message("\n▶ Processing score_min = ", sc)
    
    gene_lists_filtered <- disgenet_object[[paste0("geneLists_scoreMin", sc)]] %>%
      filter_min_vector_length(min_len = min_len)
    
    chi2_res <- run_full_overlap_analysis(
      gene_lists = c(flat_signatures[sig_names], gene_lists_filtered),
      total_genes = total_genes,
      cols_to_filter = sig_names,
      rows_to_filter = names(gene_lists_filtered),
      plot_title_or = paste0("Score: ", sc),
      triangle_mode = "full",
      fdr_threshold = 1,
      data_type = "original_data",
      verbose = FALSE
    )
    
    raw_chi2_results[[as.character(sc)]] <- chi2_res
    
    summary_by_score_chi2[[as.character(sc)]] <-
      chi2_res$processed$original_data$df %>%
      summarise_overlap() %>%
      mutate(score_min = sc)
  }
  
  summary_long_chi2 <- bind_rows(summary_by_score_chi2)
  
  # ---------------- tables ----------------
  t1_genes_all <- summary_long_chi2 %>%
    select(score_min, Var2, n_genes_all) %>%
    pivot_wider(names_from = Var2, values_from = n_genes_all) %>%
    select(score_min, any_of(signature_order))
  
  t2_genes_p005 <- summary_long_chi2 %>%
    select(score_min, Var2, n_genes_p005) %>%
    pivot_wider(names_from = Var2, values_from = n_genes_p005) %>%
    select(score_min, any_of(signature_order))
  
  t3_genes_p001 <- summary_long_chi2 %>%
    select(score_min, Var2, n_genes_p001) %>%
    pivot_wider(names_from = Var2, values_from = n_genes_p001) %>%
    select(score_min, any_of(signature_order))
  
  t4_assoc_p005 <- summary_long_chi2 %>%
    select(score_min, Var2, n_assoc_p005) %>%
    pivot_wider(names_from = Var2, values_from = n_assoc_p005) %>%
    select(score_min, any_of(signature_order))
  
  t5_assoc_p001 <- summary_long_chi2 %>%
    select(score_min, Var2, n_assoc_p001) %>%
    pivot_wider(names_from = Var2, values_from = n_assoc_p001) %>%
    select(score_min, any_of(signature_order))
  
  return(list(
    raw_chi2_results = raw_chi2_results,
    summary_long_chi2 = summary_long_chi2,
    t1_genes_all = t1_genes_all,
    t2_genes_p005 = t2_genes_p005,
    t3_genes_p001 = t3_genes_p001,
    t4_assoc_p005 = t4_assoc_p005,
    t5_assoc_p001 = t5_assoc_p001
  ))
}


##############################################################################
# ---- analysis ----
# ##############################################################################


disgenet_score_min <- c(0.9, 0.8, 0.7, 0.6, 0.5, 0)


map_df(disgenet_score_min, ~{
  
  x <- disgenet_mentalDisorders[[paste0("geneLists_scoreMin", .x)]] %>%
    lapply(length) %>% 
    unlist(use.names = FALSE)
  
  psych::describe(x) %>% 
    as.data.frame() %>%
    mutate(score_min = .x) %>% 
    select(score_min, everything()) %>% 
    rownames_to_column() %>% 
    dplyr::select(-rowname)
})

map_df(disgenet_score_min, ~{
  
  x <- disgenet_mentalDisorders[[paste0("geneLists_scoreMin", .x)]] %>%
    filter_min_vector_length(min_len = 1) %>% 
    lapply(length) %>% 
    unlist(use.names = FALSE)
  
  psych::describe(x) %>% 
    as.data.frame() %>%
    mutate(score_min = .x) %>% 
    select(score_min, everything()) %>% 
    rownames_to_column() %>% 
    dplyr::select(-rowname) %>% 
    select(-c(vars, trimmed, mad, range, skew, kurtosis, se))
})

map_df(disgenet_score_min, ~{
  
  x <- disgenet_mentalDisorders[[paste0("geneLists_scoreMin", .x)]] %>%
    filter_min_vector_length(min_len = 5) %>% 
    lapply(length) %>% 
    unlist(use.names = FALSE)
  
  psych::describe(x) %>% 
    as.data.frame() %>%
    mutate(score_min = .x) %>% 
    select(score_min, everything()) %>% 
    rownames_to_column() %>% 
    dplyr::select(-rowname) %>% 
    select(-c(vars, trimmed, mad, range, skew, kurtosis, se))
})

map_df(disgenet_score_min, ~{
  
  x <- disgenet_mentalDisorders[[paste0("geneLists_scoreMin", .x)]] %>%
    filter_min_vector_length(min_len = 10) %>% 
    lapply(length) %>% 
    unlist(use.names = FALSE)
  
  psych::describe(x) %>% 
    as.data.frame() %>%
    mutate(score_min = .x) %>% 
    select(score_min, everything()) %>% 
    rownames_to_column() %>% 
    dplyr::select(-rowname) %>% 
    select(-c(vars, trimmed, mad, range, skew, kurtosis, se))
})


# ##############################################################################
# ---- chi2 ----
# ##############################################################################
disgenetMentalHealth_GrSignatures_overlapChi2 <- run_full_overlap_analysis(
  gene_lists = c(flat_allGrSignatures_31.10.2025[sig_names], 
                 disgenet_mentalDisorders$geneLists_scoreMin0.5 %>%
                   filter_min_vector_length(min_len = 10)
  ),
  total_genes = hgnc_symbols_vector_v110,
  cols_to_filter = sig_names,
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
  col_order_original = c("globalUp", "globalDown", "NeuralCellsUp", "NeuralCellsDown",
                         "BloodCellsUp", "BloodCellsDown", "LungCellsUp", "LungCellsDown"),
  # cols_to_filter = disgenetMentalHealth_grSignatures_vector_p0.05,
  
  
  # 🎨 skala kolorów
  color_scale_range = c(-5, 5),
  text_contrast_range = c(-30, 4),
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
  # save_to_svg = "results_v2/overlap/disgenet_overlap/figures/heatmap_AllGrSignatures_disgenetMentalHealthMinScore0_log2OR_16.01.2026.svg",
  svg_width = 10.1,
  svg_height = 9.3,
  row_names_width = unit(90, "mm"),
  col_names_height = unit(50, "mm"),
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
  

# summary results
disgenetMentalHealth_GrSignatures_overlapChi2$processed$original_data$df %>% 
  filter(gene_overlap_count > 2) %>% 
  filter(p_value < 0.05) %>% 
  select(-c(._key, sig_1, fdr, fdr_value)) %>% 
  .$overlap_genes %>% strsplit(",") %>% 
  unlist %>% 
  unique() 


pgcGrSignatures_overlapChi2$processed$original_data$df %>% 
  filter(gene_overlap_count > 2) %>% 
  filter(p_value < 0.05) %>% .$overlap_genes %>% 
  strsplit(",") %>% 
  unlist %>% unique


gene_list <- disgenet_mentalDisorders$geneLists_scoreMin0.5 %>%
  filter_min_vector_length(min_len = 10)
# calculate permutation FDR
save(
  hgnc_symbols_vector_v110,
  flat_allGrSignatures_31.10.2025,
  sig_names,
  gene_list,
  file = "data/run_RData_inputs/permutation_fdr/rdata2permuationFDR_DisGeNETMentalHealth.RData"
)


# =========================
# FILTERED
# =========================
df_filtered <- disgenetMentalHealth_GrSignatures_overlapChi2$processed$original_data$df %>% 
  filter(p_value < 0.05) %>% 
  filter(gene_overlap_count >= 3) %>% 
  filter(log2_odds_ratio > 0) %>% 
  select(-c( fdr)) %>% 
  mutate(Var2 = str_replace_all(Var2, "minusGlobalUpDown5TissuesDerivedCells_", "")) %>% 
  mutate(Var2 = str_replace_all(Var2, "global_GR_genes_globalDown5TissuesDerivedCells", "systemicDown")) %>% 
  mutate(Var2 = str_replace_all(Var2, "global_GR_genes_globalUp5TissuesDerivedCells", "systemicUp"))

# =========================
# UNFILTERED
# =========================
df_unfiltered <- disgenetMentalHealth_GrSignatures_overlapChi2$processed$original_data$df %>% 
  select(-c(fdr)) %>% 
  mutate(Var2 = str_replace_all(Var2, "minusGlobalUpDown5TissuesDerivedCells_", "")) %>% 
  mutate(Var2 = str_replace_all(Var2, "global_GR_genes_globalDown5TissuesDerivedCells", "systemicDown")) %>% 
  mutate(Var2 = str_replace_all(Var2, "global_GR_genes_globalUp5TissuesDerivedCells", "systemicUp"))

# =========================
# WRITE XLSX
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
  file = "results_v2/overlap/disgenet_overlap/tables/disgenetMentalHealth_GrSignatures_overlapChi2_15.01.2026.xlsx",
  overwrite = TRUE
)
