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
  gene_lists = c(flat_allGrSignatures_31.10.2025[sig_names], 
                 random_8geneLists,
                 gene_list),
  total_genes = hgnc_symbols_vector_v110,
  cols_to_filter = c(sig_names, names(random_8geneLists)),
  rows_to_filter = c(names(gene_list)),
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

# ---- signifGrSignatures ----
heatmap_overlap_log2OR_complex(
  data_list = pgcGrSignatures_overlapChi2$processed,
  data_type = "original_data",
  
  # 🔺 triangle mode pełny
  triangle_mode = "full",
  
  # 🔹 filtrowanie
  rows_to_filter = pgc_phenotypes_vector_p0.05,
  cols_to_filter = pgc_grSignatures_vector_p0.05,
  
  # 🎨 skala kolorów
  color_scale_range = c(-3, 3),
  text_contrast_range = c(-30, 4.9),
  palette = c("#c6d3e3", "white", "darkred"),
  
  # 🔲 prostokąty istotności
  color_rects = c("#97C426", "#2F4603"),
  
  # 📊 klastrowanie
  cluster_rows = T,
  cluster_cols = T,# ============================================================
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
  # save_to_svg = "results_v2/overlap/pgc_overlap/figures/heatmap_AllGrSignaturesPGC_log2OR.svg",
  svg_width = 10.5, 
  svg_height = 9,
  row_names_width = unit(90, "mm"),
  col_names_height = unit(50, "mm"),
  force_create_directory = TRUE
)
