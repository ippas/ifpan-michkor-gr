url <- "https://www.ebi.ac.uk/gwas/api/v2/efotraits/EFO_0004247/associations/download?includeBgTraits=false&includeChildTraits=true"
# efo <- "EFO_0006788"
efo <- "EFO_0000677"
# efo <- "EFO_0005774"
# efo <- "EFO_0005774"
# EFO_0000677
url_assoc <- sprintf(
  "https://www.ebi.ac.uk/gwas/api/v2/efotraits/%s/associations/download?includeBgTraits=false&includeChildTraits=true",
  efo
)

moodDisorders_GWASCatalog <- read_tsv(url_assoc, show_col_types = FALSE, progress = FALSE)

mapperPhenotyeps2icd10 <- read_excel_sheets("data/mapperPhenotypesICD10/manualMapping_phenotype2ICD10categories_24.01.2026.xlsx") %>% 
  set_names(c( "GWASCatalog", "DisGeNET","genebass"))

mapperPhenotyeps2icd10$GWASCatalog %>% 
  mutate(
    mapping2ICD10 = mapping2ICD10 %>%
      str_replace_all("[ \\+\\(\\);]", "_") %>%  # spacja, +, (, ), ;
      str_replace_all("_+", "_") %>%             # wiele _ → jedno _
      str_replace_all("^_|_$", "")               # usuń _ na początku/końcu
  ) %>% 
  mutate(names = paste0(mapping2ICD10, "_", phenotype)) %>% 
  select(c(phenotype, names)) -> mapperPhenotypesGWASCatalog

moodDisorders_GWASCatalog %>%
  # filter(pValue < 0.000001) %>% 
  select(mappedGenes, efoTraits, pValue) %>% 
  separate_rows(mappedGenes, sep = ",") %>% 
  filter(mappedGenes %in% hgnc_symbols_vector_v110) %>% 
  left_join(., mapperPhenotypesGWASCatalog, by = c("efoTraits" = "phenotype")) %>% 
  group_by(names) %>% 
  nest %>% 
  mutate(n_genes = map(data, ~ .x %>% .$mappedGenes %>% unique %>%  length)) %>% 
  unnest(n_genes) %>%
  filter(n_genes >= 10) %>% 
  unnest() %>% 
  select(c(names, mappedGenes)) %>% 
  { split(.$mappedGenes, .$names) } %>% 
  lapply(unique) -> gene_list


random_8geneLists <- generate_random_gene_lists_by_lengths(gene_list = hgnc_symbols_vector_v110,
                                                           lengths = c(124, 208, 147, 76, 132, 86, 208, 119))


GWASCatalogMoodDisorders_GrSignatures_overlapChi2 <- run_full_overlap_analysis(
  gene_lists = c(flat_allGrSignatures_31.10.2025[sig_names][c("global_GR_genes_globalUp5TissuesDerivedCells", "global_GR_genes_globalDown5TissuesDerivedCells")], 
                 random_8geneLists,
                 gene_list
  ),
  total_genes = hgnc_symbols_vector_v110,
  cols_to_filter = c("global_GR_genes_globalUp5TissuesDerivedCells", "global_GR_genes_globalDown5TissuesDerivedCells"),
  rows_to_filter = c(names(gene_list)),
  plot_title_or = "",
  triangle_mode = "full",
  fdr_threshold = 1,
  data_type = "original_data",
  verbose = FALSE,
  palette_or = c("#c6d3e3", "white", "darkred"),
  text_contrast_range_or = c(-30, 4.9),
)

GWASCatalogMoodDisorders_GrSignatures_overlapChi2$processed$original_data$df %>% 
  filter(gene_overlap_count > 2) %>% 
  filter(p_value < 0.05)

GWASCatalogMoodDisorders_GrSignatures_overlapChi2$processed$original_data$df %>% 
  # select(-c(sig_1, ._key)) %>% 
  filter((p_value < 0.05 & gene_overlap_count >= 3 & log2_odds_ratio > 0) | grepl("F2x_schi|F2x_psy", Var1)) %>% 
  .$Var1 %>% as.character() %>% unique() -> MoodDisorders_phenotypes_vector_p0.05

GWASCatalogMoodDisorders_GrSignatures_overlapChi2$processed$original_data$df %>% 
  # select(-c(sig_1, ._key)) %>% 
  filter((p_value < 0.05 & gene_overlap_count >= 3 & log2_odds_ratio > 0)) %>% 
  # filter(gene_overlap_count >= 3) %>% 
  # filter(log2_odds_ratio > 0) %>% 
  select(-fdr) %>% 
  mutate(Var2 = str_replace_all(Var2, "minusGlobalUpDown5TissuesDerivedCells_", ""))


heatmap_overlap_log2OR_complex(
  data_list = GWASCatalogMoodDisorders_GrSignatures_overlapChi2 $processed,
  data_type = "original_data",
  
  # 🔺 triangle mode pełny
  triangle_mode = "full",
  
  # 🔹 filtrowanie
  rows_to_filter = MoodDisorders_phenotypes_vector_p0.05,
  # cols_to_filter = pgc_grSignatures_vector_p0.05,
  
  row_order_original = c(
    "G3x_Alzheimer disease",
    "G3x_Alzheimer disease,polygenic risk score",
    "G3x_Alzheimer disease,family history of Alzheimer’s disease",
    "G3x_F0x_Alzheimer disease,dementia,family history of Alzheimer’s disease",
    "G3x_F0x_late-onset Alzheimers disease",
    "F0x_dementia",
    "F3x_F2x_bipolar disorder,schizophrenia,sex interaction measurement,major depressive disorder",
    "F2x_schizophrenia",
    "F2x_schizophrenia,intelligence,self reported educational attainment",
    "F2x_schizophrenia,sex interaction measurement",
    "F2x_psychosis"
  ),
  
  p_thresholds = c(0.05, 0.01),
  
  # 🎨 skala kolorów
  color_scale_range = c(-5, 5),
  text_contrast_range = c(-30, 4.9),
  palette = c("#c6d3e3", "white", "darkred"),
  
  # 🔲 prostokąty istotnościhttp://localhost:8600/graphics/plot_zoom_png?width=1029&height=2556
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
    "global_GR_genes_globalDown5TissuesDerivedCells" =  "systemicDown",
    "global_GR_genes_globalUp5TissuesDerivedCells" =  "systemicUp"
  ),
  row_dend_height = unit(20, "mm"),
  col_dend_height = unit(20, "mm"),
  tile_gap = 1,
  save_to_svg = "results_v2/overlap/GWASCatalog_overlap/figures/heatmap_GrSignaturesGWASCatalog_log2OR_EFO0000677_25.01.2026.svg",
  svg_width = 6.63, 
  svg_height = 7.45,
  row_names_width = unit(90, "mm"),
  col_names_height = unit(50, "mm"),
  force_create_directory = TRUE,
  col_order_original = c("systemicUp", "systemicDown"
                         # "NeuralCellsUp", "NeuralCellsDown",
                         # "BloodCellsUp", "BloodCellsDown", "LungCellsUp", "LungCellsDown"
                         ),
  
)

GWASCatalogMoodDisorders_GrSignatures_overlapChi2$processed$original_data$df %>% 
  # head %>% 
  filter(gene_overlap_count > 2) %>% 
  filter(p_value < 0.05) %>% 
  filter(odds_ratio > 0) %>% 
  filter(grepl("Neural", Var2)) %>% 
  .$overlap_genes %>% 
  strsplit(",") %>% 
  unlist %>% 
  cat(sep = "\n")

# save to file
MoodDisorders_GrSignatures_overlapChi2$processed$original_data$df %>% 
  select(-c(sig_1, ._key)) %>% 
  filter(p_value < 0.05) %>% 
  filter(gene_overlap_count  >= 3) %>% 
  filter(log2_odds_ratio > 0) %>% 
  select(-fdr) %>% 
  mutate(Var2 = str_replace_all(Var2, "minusGlobalUpDown5TissuesDerivedCells_", "")) %>% 
  mutate(Var2 = str_replace_all(Var2, "global_GR_genes_globalDown5TissuesDerivedCells", "systemicDown")) %>% 
  mutate(Var2 = str_replace_all(Var2, "global_GR_genes_globalUp5TissuesDerivedCells", "systemicUp")) %>% 
  rename(phenotypeGWASCatalog = "Var1") %>% 
  rename(grSignature = "Var2") %>% 
  write.xlsx(
    file = "results_v2/overlap/GWASCatalog_overlap/GWASCatalogMoodDisorders_GRsignatures_overlap_EFO0000677.xlsx",
    rowNames = FALSE
  )


# =========================
# Sheet 1: FILTERED (jak teraz)
# =========================
df_filtered <- GWASCatalogMoodDisorders_GrSignatures_overlapChi2$processed$original_data$df %>% 
  select(-c( fdr_value, fdr)) %>%
  filter(p_value < 0.05) %>% 
  filter(gene_overlap_count >= 3) %>% 
  filter(log2_odds_ratio > 0) %>% 
  mutate(
    Var2 = str_replace_all(Var2, "minusGlobalUpDown5TissuesDerivedCells_", ""),
    Var2 = str_replace_all(Var2, "global_GR_genes_globalDown5TissuesDerivedCells", "systemicDown"),
    Var2 = str_replace_all(Var2, "global_GR_genes_globalUp5TissuesDerivedCells", "systemicUp")
  ) %>% 
  rename(
    phenotypeGWASCatalog = Var1,
    grSignature = Var2
  )

# =========================
# Sheet 2: UNFILTERED (raw, tylko kosmetyka nazw)
# =========================
df_unfiltered <- GWASCatalogMoodDisorders_GrSignatures_overlapChi2$processed$original_data$df %>% 
  select(-c( fdr_value, fdr)) %>%
  mutate(
    Var2 = str_replace_all(Var2, "minusGlobalUpDown5TissuesDerivedCells_", ""),
    Var2 = str_replace_all(Var2, "global_GR_genes_globalDown5TissuesDerivedCells", "systemicDown"),
    Var2 = str_replace_all(Var2, "global_GR_genes_globalUp5TissuesDerivedCells", "systemicUp")
  ) %>% 
  rename(
    phenotypeGWASCatalog = Var1,
    grSignature = Var2
  )

# =========================
# Create workbook
# =========================
wb <- createWorkbook()

# Sheet 1 — filtered (main results)
addWorksheet(wb, "Filtered_results")
writeData(wb, "Filtered_results", df_filtered)
freezePane(wb, "Filtered_results", firstRow = TRUE)

# Sheet 2 — unfiltered (reference / supplementary)
addWorksheet(wb, "Unfiltered_reference")
writeData(wb, "Unfiltered_reference", df_unfiltered)
freezePane(wb, "Unfiltered_reference", firstRow = TRUE)

# Save
saveWorkbook(
  wb,
  file = "results_v2/overlap/GWASCatalog_overlap/GWASCatalogMoodDisorders_GRsystemic_overlap_EFO0000677_25.01.2026.xlsx",
  overwrite = TRUE
)
# # ##############################################################################
# # calculate permutation FDR
# # ##############################################################################
# save(
#   hgnc_symbols_vector_v110,
#   flat_allGrSignatures_31.10.2025,
#   sig_names,
#   gene_list,
#   file = "data/run_RData_inputs/permutation_fdr/rdata2permuationFDR_GWASCatalogMentalHealth.RData"
# )



