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

# mapping phenotypes to ICD10
mapperPhenotyeps2icd10 <- read_excel_sheets("data/mapperPhenotypesICD10/manualMapping_phenotype2ICD10categories_24.01.2026.xlsx") %>% 
  set_names(c( "GWASCatalog", "DisGeNET","genebass"))

# ---- mapper GWAS Catalog phenotypes ----
mapperPhenotyeps2icd10$GWASCatalog %>% 
  mutate(
    mapping2ICD10 = mapping2ICD10 %>%
      str_replace_all("[ \\+\\(\\);]", "_") %>%  # spacja, +, (, ), ;
      str_replace_all("_+", "_") %>%             # wiele _ → jedno _
      str_replace_all("^_|_$", "")               # usuń _ na początku/końcu
  ) %>% 
  mutate(names = paste0(mapping2ICD10, "_", phenotype)) %>% 
  select(c(phenotype, names)) -> mapperPhenotypesGWASCatalog

# ---- mapper genebass phenotypes ----
mapperPhenotyeps2icd10$genebass %>% 
  mutate(phenotype = str_replace_all(phenotype, " ", "_")) %>% 
  mutate(names = paste0(subcategory, "_", phenotype)) %>% 
  mutate(names = str_remove(names, "NA_")) %>% 
  mutate(names = paste0(mapping2ICD10, "_", names)) %>% 
  select(c(phenotype, names)) -> mapperPhenotypesGenebass

# ---- mapper DisGeNET phenotypes ----
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


# ##############################################################################
# ---- GWAS Catalog ----
# ##############################################################################
efo <- "EFO_0000677"
efo <- "HP_0100543"
url_assoc <- sprintf(
  "https://www.ebi.ac.uk/gwas/api/v2/efotraits/%s/associations/download?includeBgTraits=false&includeChildTraits=true",
  efo
)

moodDisorders_GWASCatalog <- read_tsv(url_assoc, show_col_types = FALSE, progress = FALSE)

moodDisorders_GWASCatalog %>%
  select(mappedGenes, efoTraits, pValue) %>% 
  separate_rows(mappedGenes, sep = ",") %>% 
  filter(mappedGenes %in% hgnc_symbols_vector_v110) %>% 
  # left_join(., mapperPhenotypesGWASCatalog, by = c("efoTraits" = "phenotype")) %>% 
  mutate(names = efoTraits) %>% 
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


gene_list_disgenet <- disgenet_mentalDisorders$geneLists_scoreMin0.5 %>%
  filter_min_vector_length(min_len = 10)

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
  filter(grepl("F9x", Var1)) %>% 
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
  slice_min(order_by = p_value, n = 1, with_ties = FALSE) %>%
  arrange(source, p_value) %>% 
  .$Var1 %>% 
  unique -> GWASCatalogDisGeNETgenebass_phenotypes



# ##############################################################################
# ---- heatmap ----
# ##############################################################################

heatmap_overlap_log2OR_complex(
  data_list = GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed,
  data_type = "original_data",
  
  # 🔺 triangle mode pełny
  triangle_mode = "full",
  
  # 🔹 filtrowanie
  rows_to_filter = GWASCatalogDisGeNETgenebass_phenotypes,
  # cols_to_filter = GWASCatalogDisGeNETgenebass_phenotypes,
  
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
