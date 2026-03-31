# diseaseClasses_MSH
disgenet_metadataDiseases$diseaseClasses_MSH %>% 
  unlist() %>% 
  as.data.frame() %>% 
  set_colnames("disease") %>% 
  mutate(
    disease_name = str_trim(str_remove(disease, "\\s*\\(.*\\)$")),
    msh_code     = str_extract(disease, "(?<=\\().*(?=\\))")
  ) %>% 
  # head() %>% 
  filter(grepl("F", msh_code)) %>% unique


disgenet_metadataDiseases %>% 
  filter(diseaseClasses_MSH == "Psychological Phenomena (F02)")

disgenet_metadataDiseases %>%
  filter(grepl("diabetes", disease_name, ignore.case = T)) %>% 
  .$diseaseClasses_MSH %>% 
  unlist %>% unique()

disgenet_metadataDiseases$diseaseClasses_MSH %>% unlist %>% unique()

disgenet_metadataDiseases %>%
  dplyr::filter(sapply(diseaseClasses_MSH, function(x) "Respiratory Tract Diseases (C08)" %in% x)) %>% 
  .$disease_name %>% 
  unique %>% 
  cat(sep = "\n")
  
disgenet_metadataDiseases %>%
  dplyr::filter(sapply(diseaseClasses_UMLS_ST, function(x) "Mental or Behavioral Dysfunction (T048)" %in% x))

disgenet_metadataDiseases$diseaseClasses_MSH


#diseaseClasses_MSH
#Behavior and Behavior Mechanisms (F01)
#Nervous System Diseases (C10)
#Mental Disorders (F03)
#Psychological Phenomena (F02)


disgenet_metadataDiseases$diseaseClasses_HPO %>% 
  unlist() %>% 
  as.data.frame() %>% 
  set_colnames("disease") %>% unique


intersect(
  disgenet_RespiratoryTractDisease$geneLists_scoreMin0$Asthma,
  disgenet_mentalDisorders$geneLists_scoreMin0$Major_Depressive_Disorder
)

disgenet_RespiratoryTractDisease$geneLists_scoreMin0 %>% 
  filter_min_vector_length(min_len = 50) %>% length()

intersect(
  disgenet_RespiratoryTractDisease$geneLists_scoreMin0$Chronic_Obstructive_Airway_Disease,
  disgenet_mentalDisorders$geneLists_scoreMin0$Major_Depressive_Disorder
)


disgenet_mentalDisorders$geneLists_scoreMin0$`Alzheimer's_Disease`


disgenet_mentalDisorders$geneLists_scoreMin0[c(
  "Major_Depressive_Disorde",
  "Schizophrenia",
  "Alcohol_abuse",
  "Bipolar_Disorder",
  "Alzheimer's_Disease"
)]

disgenet_mentalDisorders$geneLists_scoreMin0$Anxiety_Disorders


mental_sets <- disgenet_mentalDisorders$geneLists_scoreMin0[c(
  "Major_Depressive_Disorder",
  "Unipolar_Depression",
  "Depressive_disorder",
  "Bipolar_Disorder",
  "Schizophrenia",
  "Alzheimer's_Disease",
  "Anxiety_Disorders",
  "Autism_Spectrum_Disorders"   # <- dopasuj do realnej nazwy u Ciebie
)]

resp_sets <- disgenet_RespiratoryTractDisease$geneLists_scoreMin0 %>%
  filter_min_vector_length(min_len = 10)

mental_sets <- mental_sets[!is.na(names(mental_sets))]
mental_sets <- mental_sets[!vapply(mental_sets, is.null, logical(1))]

mentalLungPhenotypes_overlapChi2 <- run_full_overlap_analysis(
  gene_lists      = c(mental_sets, resp_sets),
  total_genes     = hgnc_symbols_vector_v110,
  cols_to_filter  = names(mental_sets),
  rows_to_filter  = names(resp_sets),
  plot_title_or   = "",
  triangle_mode   = "full",
  fdr_threshold   = 1,
  data_type       = "original_data",
  verbose         = FALSE,
  palette_or      = c("#c6d3e3", "white", "darkred"),
  text_contrast_range_or = c(-30, 4.9)
)


heatmap_overlap_log2OR_complex(
  data_list = mentalLungPhenotypes_overlapChi2$processed,
  data_type = "original_data",
  
  # 🔺 triangle mode pełny
  triangle_mode = "full",
  
  # 🔹 filtrowanie
  # rows_to_filter = GWASCatalogDisGeNETgenebass_phenotypes,
  
  # row_order_original = GWASCatalogDisGeNETgenebass_phenotypes,
  
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
  # col_mapper = c(
  #   "global_GR_genes_globalDown5TissuesDerivedCells" =  "systemicDown",
  #   "global_GR_genes_globalUp5TissuesDerivedCells" =  "systemiclUp"
  # ),
  row_dend_height = unit(20, "mm"),
  col_dend_height = unit(20, "mm"),
  tile_gap = 1,
  # save_to_svg = "results_v2/overlap/overlap_GWASCatalogDisGeNETgenebass/figures/heatmap_GRsystemic_GWASCatalogDisGeNETgenebass_log2OR_v2_28.01.2026.svg",
  svg_width = 6.625,
  svg_height = 22,
  row_names_width = unit(90, "mm"),
  col_names_height = unit(50, "mm"),
  force_create_directory = TRUE
)

mentalLungPhenotypes_overlapChi2$processed$original_data$df %>% 
  filter(gene_overlap_count >= 20 ) %>% 
  .$overlap_genes %>% 
  strsplit(",") %>% unlist() %>% 
  unique()

mentalLungPhenotypes_overlapChi2$processed$original_data$df %>% 
  filter(Var2 == "Major_Depressive_Disorder") %>% 
  filter(grepl("CYP2D6|DTNBP1|EHD3|HSPA1B|NTRK2|PAWR", overlap_genes))

mentalLungPhenotypes_overlapChi2$processed$original_data$df %>% 
  filter(Var2 == "Major_Depressive_Disorder") %>% 
  mutate(
    overlapping_selected_genes = sapply(
      strsplit(overlap_genes, ",\\s*"),
      function(x) paste(x[grepl("CYP2D6|DTNBP1|EHD3|HSPA1B|NTRK2|PAWR", x)], collapse = ", ")
    )
  ) %>% 
  filter(overlapping_selected_genes != "")

mentalLungPhenotypes_overlapChi2$processed$original_data$df %>% 
  filter(Var2 == "Schizophrenia") %>% 
  mutate(
    overlapping_selected_genes = sapply(
      strsplit(overlap_genes, ",\\s*"),
      function(x) paste(x[grepl("ATXN1|BCL2|IL1B|JUN|L1CAM|LRP8|MYO5B|NR4A2|NRG1|PDE4B|PLA2G4A|PLAT|RTKN2|SEMA3A|SLC6A9|SOCS2", x)], collapse = ", ")
    )
  ) %>% 
  filter(overlapping_selected_genes != "")

mentalLungPhenotypes_overlapChi2$processed$original_data$df %>% 
  filter(Var2 == "Bipolar_Disorder") %>% 
  filter(gene_overlap_count >= 10 ) %>% 
  mutate(
    overlapping_selected_genes = sapply(
      strsplit(overlap_genes, ",\\s*"),
      function(x) paste(x[grepl("BCL2|BDKRB2|DIO2|IL1B|NR4A2|NRG1|PDE4B|PLA2G4A|SNAP25", x)], collapse = ", ")
    )
  ) %>% 
  filter(overlapping_selected_genes != "")



disgenet_cardiovascularDiseasesC14 <- readRDS("data/databases/disgenet/disgenet_MSH_cardiovascularDiseasesC14.rds")
disgenet_immuneSystemDiseasesC20 <- readRDS("data/databases/disgenet/disgenet_MSH_ImmuneSystemDiseasesC20.rds")
disgenet_nutritionalAndMetabolicDiseasesC18 <- readRDS("data/databases/disgenet/disgenet_MSH_NutritionalAndMetabolicDiseasesC18.rds")

disgenet_cardiovascularDiseasesC14$geneLists_scoreMin0 %>% 
  filter_min_vector_length(50) %>% length()



mentalLungPhenotypes_overlapChi2 <- run_full_overlap_analysis(
  gene_lists      = c(mental_sets,
                      flat_allGrSignatures_31.10.2025[sig_names],
                      disgenet_immuneSystemDiseasesC20$geneLists_scoreMin0 %>% 
                        filter_min_vector_length(10)),
  total_genes     = hgnc_symbols_vector_v110,
  cols_to_filter  = c(names(mental_sets), sig_names),
  rows_to_filter  = names(disgenet_immuneSystemDiseasesC20$geneLists_scoreMin0 %>% 
                            filter_min_vector_length(10)),
  plot_title_or   = "",
  triangle_mode   = "full",
  fdr_threshold   = 1,
  data_type       = "original_data",
  verbose         = FALSE,
  palette_or      = c("#c6d3e3", "white", "darkred"),
  text_contrast_range_or = c(-30, 4.9)
)



heatmap_overlap_log2OR_complex(
  data_list = mentalLungPhenotypes_overlapChi2$processed,
  data_type = "original_data",
  
  # 🔺 triangle mode pełny
  triangle_mode = "full",
  
  # 🔹 filtrowanie
  # rows_to_filter = GWASCatalogDisGeNETgenebass_phenotypes,
  
  # row_order_original = GWASCatalogDisGeNETgenebass_phenotypes,
  
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
  # col_mapper = c(
  #   "global_GR_genes_globalDown5TissuesDerivedCells" =  "systemicDown",
  #   "global_GR_genes_globalUp5TissuesDerivedCells" =  "systemiclUp"
  # ),http://localhost:8600/graphics/plot_zoom_png?width=1180&height=1440
  row_dend_height = unit(20, "mm"),
  col_dend_height = unit(20, "mm"),
  tile_gap = 1,
  # save_to_svg = "results_v2/overlap/overlap_GWASCatalogDisGeNETgenebass/figures/heatmap_GRsystemic_GWASCatalogDisGeNETgenebass_log2OR_v2_28.01.2026.svg",
  svg_width = 6.625,
  svg_height = 22,
  row_names_width = unit(90, "mm"),
  col_names_height = unit(50, "mm"),
  force_create_directory = TRUE
)
