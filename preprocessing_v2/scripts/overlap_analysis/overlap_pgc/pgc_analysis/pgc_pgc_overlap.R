# ##############################################################################
# ---- uses data ----
# ##############################################################################
pgc_geneList_10e4$pgc_PGC3_SCZ_wave3.asian.autosome.public.v3.vcf.tsv.tsv

disgenet_mentalDisorders$geneLists_scoreMin0.5$Schizophrenia


intersect(pgc_geneList_10e4$pgc_PGC3_SCZ_wave3.asian.autosome.public.v3.vcf.tsv.tsv,disgenet_mentalDisorders$geneLists_scoreMin0.5$Schizophrenia)



# ============================================================
# 🧩 2️⃣ Uruchomienie analizy overlap (dla nowych sygnatur)
# ============================================================

pgc_overlapChi2 <- run_full_overlap_analysis(
  gene_lists = pgc_geneList_10e4,
  total_genes = hgnc_symbols_vector_v110,
  cols_to_filter = names(pgc_geneList_10e4),
  rows_to_filter = names(pgc_geneList_10e4),
  plot_title_or = "",
  triangle_mode = "full",
  fdr_threshold = 1,
  data_type = "original_data",
  verbose = FALSE
)



# ##############################################################################
# ---- complex heatmap ----
# ##############################################################################

# ---- All GrSignatures ----
heatmap_overlap_log2OR_complex(
  data_list = pgc_overlapChi2$processed,
  data_type = "original_data",
  
  # 🔺 triangle mode pełny
  triangle_mode = "full",
  
  # 🔹 filtrowanie
  # rows_to_filter = pgc_phenotypes_vector_p0.05,
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
