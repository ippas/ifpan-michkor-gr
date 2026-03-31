# ============================================================
# 📚 Pakiety
# ============================================================
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)

# ============================================================
# 🧬 1️⃣ Przygotowanie listy genów PGC (p < 1e-4)
# ============================================================

gene_list <- pgc_annotation_geneCenter50kb_p1e4 %>%
  filter(pvalue < 0.0001) %>%
  select(gene_symbol, source_file) %>%
  distinct() %>%
  group_by(source_file) %>%
  summarise(genes = list(unique(gene_symbol)), .groups = "drop") %>%
  deframe()


gene_list <- pgc_annotation_geneCenter50kb_p1e4 %>%
  filter(pvalue < 0.0001) %>%
  select(gene_symbol, source_file) %>% 
  filter(gene_symbol %in% hgnc_symbols_vector_v110) %>% 
  unique %>%
  group_by(source_file) %>% 
  nest() %>% 
  mutate(n_genes = map(data, ~ .x %>% nrow)) %>% 
  unnest(n_genes) %>% 
  filter(n_genes >= 10) %>% 
  select(-n_genes) %>% 
  unnest(data) %>% 
  ungroup %>% 
  distinct() %>%
  group_by(source_file) %>%
  summarise(genes = list(unique(gene_symbol)), .groups = "drop") %>%
  deframe()

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


pgcGrSignatures_overlapChi2 <- run_full_overlap_analysis(
  gene_lists = c(flat_allGrSignatures_31.10.2025[sig_names], gene_list),
  total_genes = hgnc_symbols_vector_v110,
  cols_to_filter = sig_names,
  rows_to_filter = names(gene_list),
  plot_title_or = "",
  triangle_mode = "full",
  fdr_threshold = 1,
  data_type = "original_data",
  verbose = FALSE
)

# ============================================================
# 📊 3️⃣ Szybkie sprawdzenie rozmiaru wyników
# ============================================================

cat("Significant data dimensions: ",
    paste(dim(pgcGrSignatures_overlapChi2$processed$significant_data$df), collapse = " x "), "\n")
cat("Original data dimensions: ",
    paste(dim(pgcGrSignatures_overlapChi2$processed$original_data$df), collapse = " x "), "\n")

# ============================================================
# 🧮 4️⃣ Filtrowanie i wyciągnięcie genów z istotnych overlapów
# ============================================================

tmp <- pgcGrSignatures_overlapChi2

genes_significant <- tmp$processed$original_data$df %>%
  filter(
    Var1 %in% names(gene_list),
    Var2 %in% sig_names
  ) %>%
  select(-fdr) %>%
  group_by(Var1) %>%
  nest() %>%
  mutate(
    data = map(data, ~ .x %>%
                 mutate(fdr = p.adjust(p_value, method = "fdr")))
  ) %>%
  unnest(data) %>%
  filter(
    gene_overlap_count > 2,
    log2_odds_ratio > 0,
    fdr < 0.1
  ) %>%
  pull(overlap_genes) %>%
  strsplit(",") %>%
  unlist() %>%
  unique()

# ============================================================
# 📈 5️⃣ Liczba unikalnych genów GR-zależnych w istotnych overlapach
# ============================================================

length(genes_significant)

pgcGrSignatures_overlapChi2$processed$original_data$df %>%
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
  filter(p_value < 0.05) -> pgcGrSignatures_filtered_df

# ============================================================
# 🧠 4️⃣ Wyznacz listę nazw sygnatur z istotnymi overlapami
# ============================================================

pgcGrSignatures_filtered_df$Var2 %>% unique() -> pgc_grSignatures_vector_p0.05

pgcGrSignatures_filtered_df$Var1 %>% unique() -> pgc_phenotypes_vector_p0.05

pgcGrSignatures_filtered_df %>% 
  filter(Var1 %in% c(
    "pgc_mdd_symptoms_2023-Comm-MDD7_worthless.txt.tsv",
    "pgc_pts_eur_freeze2_overall.results.tsv",
    "pgc_mdd_symptoms_2023-Comm-MDD9_death.txt.tsv",
    "pgc_mdd_symptoms_2023-Comm-MDD1_depressed.txt.tsv",
    "pgc_PGC_MDD2018_10kSNPs.tsv",
    "pgc_mdd_symptoms_2023-Clin-MDD9_death.txt.tsv",
    "pgc_TS_Oct2018.tsv"
  )) %>% .$Var1 %>% unique() -> pgc_phenotypes_vector_p0.05
# ============================================================
# 🧮 5️⃣ Liczba genów dla każdej sygnatury przy różnych progach
# ============================================================

pgc_summary_df <- pgcGrSignatures_overlapChi2$processed$original_data$df %>%
  filter(Var2 %in% sig_names) %>%
  group_by(Var2) %>%
  nest() %>%
  mutate(
    data = map(data, ~ .x %>%
                 mutate(fdr = p.adjust(p_value, method = "fdr")))
  ) %>%
  mutate(
    n_genes_all = map_int(data, ~ .x %>%
                            filter(overlap_genes != "") %>%
                            pull(overlap_genes) %>%
                            strsplit(",") %>%
                            unlist() %>%
                            unique() %>%
                            length()),
    n_genes_p0.05 = map_int(data, ~ .x %>%
                              filter(p_value < 0.05,
                                     log2_odds_ratio > 0,
                                     overlap_genes != "",
                                     gene_overlap_count > 2) %>%
                              pull(overlap_genes) %>%
                              strsplit(",") %>%
                              unlist() %>%
                              unique() %>%
                              length()),
    n_genes_p0.01 = map_int(data, ~ .x %>%
                              filter(p_value < 0.01,
                                     log2_odds_ratio > 0,
                                     overlap_genes != "",
                                     gene_overlap_count > 2) %>%
                              pull(overlap_genes) %>%
                              strsplit(",") %>%
                              unlist() %>%
                              unique() %>%
                              length())
  ) %>%
  select(Var2, n_genes_all, n_genes_p0.05, n_genes_p0.01)

# ============================================================
# 📊 6️⃣ Wizualizacja overlapów (log2(OR))
# ============================================================

heatmap_overlap_log2OR_ggplot(
  data_list = pgcGrSignatures_overlapChi2$processed,
  data_type = "original_data",
  triangle_mode = "full",
  title = "Odds ratio overlap (PGC vs GR signatures)",
  color_scale_range = c(-3, 3),
  text_contrast_range = c(-30, 4.9),
  palette = c("#c6d3e3", "white", "darkred"),
  rows_to_filter = pgc_phenotypes_vector_p0.05,
  # cols_to_filter =  pgc_grSignatures_vector_p0.05,
  color_rects = c("#97C426", "#2F4603"),
  cluster_rows = F,
  cluster_cols = F,
  
  show_dendrograms = TRUE,
  sig_rect_linewidth = 1
) -> p1

custom_legend <- create_customRect_patch_legend(
  colors = c("#97C426", "#2F4603"),
  labels = c("p < 0.05", "p < 0.01"),
  spacing = 4,
  box_linewidth = 4
)

final_plot <- wrap_plots(p1, custom_legend, ncol = 1, heights = c(10, 1))

final_plot



# ##############################################################################
# ---- complex heatmap ----
# ##############################################################################

# ---- All GrSignatures ----
heatmap_overlap_log2OR_complex(
  data_list = pgcGrSignatures_overlapChi2$processed,
  data_type = "original_data",
  
  # 🔺 triangle mode pełny
  triangle_mode = "full",
  
  # 🔹 filtrowanie
  rows_to_filter = pgc_phenotypes_vector_p0.05,
  # cols_to_filter = pgc_grSignatures_vector_p0.05,
  
  # 🎨 skala kolorów
  color_scale_range = c(-3, 3),
  text_contrast_range = c(-30, 4.9),
  palette = c("#c6d3e3", "white", "darkred"),
  
  # 🔲 prostokąty istotności
  color_rects = c("#97C426", "#2F4603"),
  
  # 📊 klastrowanie
  cluster_rows = T,
  cluster_cols = T,
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


# ---- signifGrSignatures ----
heatmap_overlap_log2OR_complex(
  data_list = pgcGrSignatures_overlapChi2$processed,
  data_type = "original_data",
  
  # 🔺 triangle mode pełny
  triangle_mode = "full",
  
  # 🔹 filtrowanie
  rows_to_filter = pgc_phenotypes_vector_p0.05,
  # cols_to_filter = pgc_grSignatures_vector_p0.05,
  
  # 🎨 skala kolorów
  color_scale_range = c(-5, 5),
  text_contrast_range = c(-30, 4.9),
  palette = c("#c6d3e3", "white", "darkred"),
  
  # 🔲 prostokąty istotności
  color_rects = c("#97C426", "#2F4603"),
  
  # 📊 klastrowanie
  cluster_rows = T,
  # cluster_cols = T,
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
  save_to_svg = "results_v2/overlap/pgc_overlap/figures/heatmap_SignifGrSignaturesPGC_log2OR_13.01.2026.svg",
  svg_width = 10.5, 
  svg_height = 6.2,
  row_names_width = unit(90, "mm"),
  col_names_height = unit(50, "mm"),
  force_create_directory = TRUE,
  col_order_original = c("globalUp", "globalDown", "NeuralCellsUp", "NeuralCellsDown", "BloodCellsUp", "BloodCellsDown", "LungCellsUp", "LungCellsDown")
)



pgcGrSignatures_filtered_df %>% 
  filter(Var1 %in% c(
    "pgc_mdd_symptoms_2023-Comm-MDD7_worthless.txt.tsv",
    "pgc_pts_eur_freeze2_overall.results.tsv",
    "pgc_mdd_symptoms_2023-Comm-MDD9_death.txt.tsv",
    "pgc_mdd_symptoms_2023-Comm-MDD1_depressed.txt.tsv",
    "pgc_PGC_MDD2018_10kSNPs.tsv",
    "pgc_mdd_symptoms_2023-Clin-MDD9_death.txt.tsv",
    "pgc_TS_Oct2018.tsv"
  )) %>% 
  select(-c(fdr_value, fdr, sig_1, ._key)) %>% 
  mutate(Var2 = str_replace_all(Var2, "minusGlobalUpDown5TissuesDerivedCells_", "")) %>% 
  mutate(Var2 = str_replace_all(Var2, "global_GR_genes_globalDown5TissuesDerivedCells", "globalDown")) %>% 
  rename(grSignature = "Var2") %>% 
  rename(GWASfile = "Var1") %>% 
  write.xlsx("results_v2/overlap/pgc_overlap/tables//PGC_GrSignatures_overlapChi2_13.01.2026.xlsx", 
             rowNames = FALSE)

