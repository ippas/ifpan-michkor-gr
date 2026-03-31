# ============================================
# 📦 Required packages
# ============================================
if (!requireNamespace("ontologyIndex", quietly = TRUE))
  install.packages("ontologyIndex")

library(ontologyIndex)
library(dplyr)
library(stringr)

# ============================================
# 📁 Local path for ontology storage
# ============================================
local_dir <- "/home/mateusz/projects/ifpan-michkor-gr/data/databases/mgi_phenotypes"
mp_path   <- file.path(local_dir, "mp.obo")

# ============================================
# 🧬 Function: Find MP descendants by phenotype name
# ============================================
mp_get_descendants_by_name <- function(category_name, mp = NULL, file_path = mp_path) {
  
  # Load ontology if not provided
  if (is.null(mp)) {
    if (!file.exists(file_path)) {
      message("Downloading mp.obo ontology file...")
      download.file(
        url      = "http://purl.obolibrary.org/obo/mp.obo",
        destfile = file_path,
        quiet    = TRUE
      )
    } else {
      message("mp.obo already exists — using local ontology file")
    }
    mp <- get_ontology(file_path, extract_tags = "everything")
  }
  
  # Find category ID by name (case-insensitive)
  category_id <- mp$id[tolower(mp$name) == tolower(category_name)]
  if (length(category_id) == 0) {
    stop("Phenotype category not found: ", category_name)
  }
  
  # Retrieve all descendant phenotype terms
  descendants <- get_descendants(mp, category_id)
  
  # Build output dataframe
  df <- data.frame(
    term_id = descendants,
    name = mp$name[descendants],
    stringsAsFactors = FALSE
  )
  
  rownames(df) <- NULL
  message("✔ Found ", nrow(df), " descendant phenotypes for: ", category_name)
  return(df)
}

# ============================================
# 🎯 Example usage
# ============================================
behaviorNeurological_mgi <- mp_get_descendants_by_name(
  category_name = "cardiovascular system phenotype"
)

# ============================================
# 🔗 Download MGI 2024 phenotype gene associations
# ============================================
mgi_url <- "https://raw.githubusercontent.com/MaayanLab/Enrichr-Viz-Appyter/master/Enrichr-Processed-Library-Storage/Clustered_Scatterplots/MGI_Mammalian_Phenotype_Level_4_2024.csv"
mgi2024 <- read.csv(mgi_url, stringsAsFactors = FALSE)

# ============================================
# 🧪 Extract gene sets for filtered MP terms
# ============================================
behaviorNeurological_genes <- mgi2024 %>%
  select(term, genes) %>%
  mutate(mp_id = str_extract(term, "MP:\\d{7}")) %>%
  filter(mp_id %in% behaviorNeurological_mgi$term_id)

behaviorNeurological_genes %>%
  separate_rows(genes, sep = " ") %>%     # split genes by space
  filter(genes != "") %>%                 # remove empty tokens if any
  arrange(mp_id, genes) %>% 
  select(c(term, mp_id, genes)) %>% 
  filter(genes %in% hgnc_symbols_vector_v110) %>%
  group_by(term, mp_id) %>% 
  nest %>% 
  mutate(n_genes = map(data, ~ .x %>% nrow)) %>% 
  unnest(n_genes) %>% 
  filter(n_genes >= 10) %>%
  select(-n_genes) %>% 
  select(-mp_id) %>% 
  unnest(data) %>% 
  ungroup %>% 
  distinct() %>%
  group_by(term) %>% 
  summarise(genes = list(unique(genes)), .groups = "drop") %>%
  deframe() -> gene_list

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

GWASCatalog_GrSignatures_overlapChi2 <- run_full_overlap_analysis(
  gene_lists = c(flat_allGrSignatures_18.11.2025[sig_names], gene_list),
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

GWASCatalog_GrSignatures_overlapChi2$processed$original_data$df %>% head
  filter(p_value < 0.05) %>% 
  filter(gene_overlap_count  >= 3) %>% 
  filter(log2_odds_ratio > 0) %>% 
  .$Var2 %>% as.character() %>% unique() -> mgi_grSignatures_vector_p0.05

mgiGrSignatures_overlapChi2$processed$original_data$df %>% 
  filter(p_value < 0.05) %>% 
  filter(gene_overlap_count  >= 3) %>%
  filter(log2_odds_ratio > 0) %>% 
  .$Var1 %>% as.character() %>% unique() -> mgi_phenotypes_vector_p0.05

# ##############################################################################
# ---- complex heatmap ----
# ##############################################################################

# ---- All GrSignatures ----
heatmap_overlap_log2OR_complex(
  data_list = mgiGrSignatures_overlapChi2$processed,
  data_type = "original_data",
  
  # 🔺 triangle mode pełny
  triangle_mode = "full",
  
  # 🔹 filtrowanie
  rows_to_filter = mgi_phenotypes_vector_p0.05,
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





  
  