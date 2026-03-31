disgenet_mentalDisorders$geneLists_scoreMin0.5 %>%
  filter_min_vector_length(min_len = 10)



disgenet_overlapChi2 <- run_full_overlap_analysis(
  gene_lists = disgenet_mentalDisorders$geneLists_scoreMin0.5 %>%
    filter_min_vector_length(min_len = 10),
  total_genes = hgnc_symbols_vector_v110,
  cols_to_filter = names(disgenet_mentalDisorders$geneLists_scoreMin0.5 %>%
                           filter_min_vector_length(min_len = 10)),
  rows_to_filter = names(disgenet_mentalDisorders$geneLists_scoreMin0.5 %>%
                           filter_min_vector_length(min_len = 10)),
  plot_title_or = "",
  triangle_mode = "full",
  fdr_threshold = 1,
  data_type = "original_data",
  verbose = FALSE
)


# ============================================================
# 🧠 4️⃣ Wyznacz listę nazw sygnatur z istotnymi overlapami
# ============================================================

disgenet_overlapChi2$processed$original_data$df %>% filter(p_value < 0.05) %>% 
  filter(gene_overlap_count >= 10) %>% .$Var1 %>% unique -> disgenet_vector_p0.05


# ##############################################################################
# ---- complex heatmap ----
# ##############################################################################

# ---- All GrSignatures ----
heatmap_overlap_log2OR_complex(
  data_list = disgenet_overlapChi2$processed,
  data_type = "original_data",
  
  # 🔺 triangle mode pełny
  triangle_mode = "full",
  
  # 🔹 filtrowanie
  rows_to_filter = disgenet_vector_p0.05,
  cols_to_filter = disgenet_vector_p0.05,
  
  # 🎨 skala kolorów
  color_scale_range = c(-8, 8),
  text_contrast_range = c(-30, 4.9),
  palette = c("#c6d3e3", "white", "darkred"),
  
  # 🔲 prostokąty istotności
  color_rects = c("#97C426", "#2F4603"),
  
  # 📊 klastrowanie
  cluster_rows = T,
  cluster_cols = F,
  show_dendrograms = TRUE,
  rect_lwd = 2.5,

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

