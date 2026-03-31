
gene_list_disgenet
gene_list_genebass
gene_list_GWASCatalog


gtex_tmp <- multi_downloadMedianSubtissuesGTEx_tissueSummaryExpression(gene_symbol = gene_list_disgenet$DisGeNET_F0x_Amnesia,
                                                                       verbose = F,
                                                                       show_progress = T,
                                                                       quiet_single = F)

gtex_tmp %>% 
  filter(grepl("Brain|Lung|Whole_Blood", tissue)) %>% 
  filter(median_max > 1) -> gtex_tmp



run_GTEx_batch <- function(gene_list_object,
                           tissue = "Whole_Blood",
                           median_max_threshold = 1) {
  
  t_start <- Sys.time()
  
  # 1) Unia genów (unikalne)
  all_genes <- unique(unlist(gene_list_object, use.names = FALSE))
  all_genes <- all_genes[!is.na(all_genes) & all_genes != ""]
  
  if (length(all_genes) == 0) {
    stop("No genes to query: union of gene_list_object is empty.")
  }
  
  cat(sprintf("Unique genes to query: %d\n", length(all_genes)))
  cat(sprintf("Downloading GTEx once... (target tissue: %s)\n", tissue))
  
  # 2) Pobierz GTEx raz
  gtex_all <- multi_downloadMedianSubtissuesGTEx_tissueSummaryExpression(
    gene_symbol = all_genes
  )
  
  # 3) Geny przechodzące filtr dla wybranej tkanki
  # (u Ciebie kolumna z genami to geneSymbol)
  if (!("geneSymbol" %in% names(gtex_all))) {
    stop(
      "Column 'geneSymbol' not found in GTEx output. Available columns: ",
      paste(names(gtex_all), collapse = ", ")
    )
  }
  if (!("tissue" %in% names(gtex_all))) {
    stop(
      "Column 'tissue' not found in GTEx output. Available columns: ",
      paste(names(gtex_all), collapse = ", ")
    )
  }
  if (!("median_max" %in% names(gtex_all))) {
    stop(
      "Column 'median_max' not found in GTEx output. Available columns: ",
      paste(names(gtex_all), collapse = ", ")
    )
  }
  
  gtex_pass_genes <- gtex_all %>%
    dplyr::filter(tissue == tissue) %>%                     # celowo: nazwa argumentu i kolumny różne? patrz niżej fix
    dplyr::filter(median_max > median_max_threshold) %>%
    dplyr::pull(geneSymbol) %>%
    unique()
  
  # IMPORTANT FIX: unikamy konfliktu nazwy 'tissue' (argument) vs 'tissue' (kolumna)
  # Jeśli chcesz mieć to od razu poprawnie, podmień powyższy filter na:
  # dplyr::filter(.data$tissue == tissue)
  
  # Poprawiona wersja filtra:
  gtex_pass_genes <- gtex_all %>%
    dplyr::filter(.data$tissue == tissue) %>%
    dplyr::filter(.data$median_max > median_max_threshold) %>%
    dplyr::pull(.data$geneSymbol) %>%
    unique()
  
  cat(sprintf(
    "Genes passing filter in %s (median_max > %s): %d\n",
    tissue, median_max_threshold, length(gtex_pass_genes)
  ))
  
  # 4) Przefiltruj wejściowy obiekt -> lista wektorów genów
  out <- vector("list", length(gene_list_object))
  names(out) <- names(gene_list_object)
  
  n_total <- length(gene_list_object)
  i <- 0
  
  for (nm in names(gene_list_object)) {
    i <- i + 1
    cat(sprintf("[%d/%d | %.1f%%] %s\n", i, n_total, 100 * i / n_total, nm))
    
    genes <- gene_list_object[[nm]]
    genes <- genes[!is.na(genes) & genes != ""]
    
    out[[nm]] <- intersect(genes, gtex_pass_genes)
    cat(sprintf("  ✓ kept: %d / %d\n", length(out[[nm]]), length(genes)))
  }
  
  # 5) Czas wykonania
  elapsed_sec <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  cat(sprintf("\nDone. Total runtime: %.1f seconds (%.2f minutes)\n",
              elapsed_sec, elapsed_sec / 60))
  
  return(out)
}

disgenet_genes_expressed_brain <- run_GTEx_batch(
  gene_list_disgenet,
  tissue = "Brain",
  median_max_threshold = 1
)


GWASCatalog_genes_expressed_brain <- run_GTEx_batch(
  gene_list_GWASCatalog,
  tissue = "Brain",
  median_max_threshold = 5
)

gene_list_GWASCatalog %>% unname %>% unlist %>% unique() %>% length()

GWASCatalog_genes_expressed_brain %>% unname %>% unlist %>% unique() %>% length()

# ---- test overlap ----

sig_names <- c(
  # "minusGlobalUpDown5TissuesDerivedCells_BloodCellsUp",
  # "minusGlobalUpDown5TissuesDerivedCells_LungCellsUp",
  "minusGlobalUpDown5TissuesDerivedCells_NeuralCellsUp",
  # "minusGlobalUpDown5TissuesDerivedCells_BloodCellsDown",
  # "minusGlobalUpDown5TissuesDerivedCells_LungCellsDown",
  "minusGlobalUpDown5TissuesDerivedCells_NeuralCellsDown"
)


sandobx_overlap <- run_full_overlap_analysis(
  gene_lists = c(flat_allGrSignatures_31.10.2025[sig_names], 
                 GWASCatalog_genes_expressed_brain
  ),
  total_genes = hgnc_symbols_vector_v110,
  cols_to_filter = c(sig_names),
  rows_to_filter = c(names(GWASCatalog_genes_expressed_brain)),
  plot_title_or = "",
  triangle_mode = "full",
  fdr_threshold = 1,
  data_type = "original_data",
  verbose = FALSE,
  palette_or = c("#c6d3e3", "white", "darkred"),
  text_contrast_range_or = c(-30, 4.9)
)





sandobx_overlap$processed$original_data$df %>% 
  # filter(grepl("F1x|F2x|F3x|F4x", Var1)) %>% 
  # filter(grepl("F4x", Var1)) %>% 
  filter((p_value < 0.05 & gene_overlap_count >= 3 & log2_odds_ratio > 0)) %>%
  # filter(gene_overlap_count >= 3) %>%
  # filter(log2_odds_ratio > 0) %>%
  .$Var1 %>% as.character() %>% unique()  -> genebass_phenotypes_vector_p0.05

# ##############################################################################
# ---- complex heatmap ----
# ##############################################################################

# ---- All GrSignatures ----
heatmap_overlap_log2OR_complex(
  data_list = sandobx_overlap$processed,
  data_type = "original_data",
  
  # 🔺 triangle mode pełny
  triangle_mode = "full",
  
  # 🔹 filtrowanie
  rows_to_filter = genebass_phenotypes_vector_p0.05,
  
  col_order_original = c(
    # "globalUp", "globalDown" 
                         "NeuralCellsUp", "NeuralCellsDown"
                         # "BloodCellsUp", "BloodCellsDown", "LungCellsUp", "LungCellsDown"
  ),

  # cols_to_filter = pgc_grSignatures_vector_p0.05,
  
  # 🎨 skala kolorów
  color_scale_range = c(-5, 5),
  text_contrast_range = c(-30, 4.9),
  palette = c("#c6d3e3", "white", "darkred"),
  
  # 🔲 prostokąty istotności
  color_rects = c("#97C426", "#2F4603"),
  
  # 📊 klastrowanie
  cluster_rows = F,
  cluster_cols = F,
  show_dendrograms = TRUE,
  rect_lwd = 2.5,
  col_mapper = c(
    # "minusGlobalUpDown5TissuesDerivedCells_BloodCellsUp" = "BloodCellsUp",
    # "minusGlobalUpDown5TissuesDerivedCells_LungCellsUp" = "LungCellsUp",
    "minusGlobalUpDown5TissuesDerivedCells_NeuralCellsUp" = "NeuralCellsUp",
    # "minusGlobalUpDown5TissuesDerivedCells_BloodCellsDown" = "BloodCellsDown",
    # "minusGlobalUpDown5TissuesDerivedCells_LungCellsDown" = "LungCellsDown",
    "minusGlobalUpDown5TissuesDerivedCells_NeuralCellsDown" = "NeuralCellsDown"
    # "global_GR_genes_globalDown5TissuesDerivedCells" =  "globalDown",
    # "global_GR_genes_globalUp5TissuesDerivedCells" =  "globalUp"
  ),
  row_dend_height = unit(20, "mm"),
  col_dend_height = unit(20, "mm"),
  tile_gap = 1,
  # save_to_svg = "results_v2/overlap/genebass_overlap/figures/heatmap_systemicGrSignaturesGenebass_log2OR_25.01.2026.svg",
  svg_width = 6.62, 
  svg_height = 7.25,
  row_names_width = unit(90, "mm"),
  col_names_height = unit(50, "mm"),
  force_create_directory = TRUE
)


GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$df %>% 
  filter(p_value < 0.05, gene_overlap_count >= 3, log2_odds_ratio > 0) %>% 
  select(-c(fdr, fdr_value)) %>% 
  filter(grepl("_F3x_", Var1)) %>% 
  filter(!grepl("F2x", Var1)) %>% 
  .$overlap_genes %>% 
  strsplit(",") %>% 
  unlist %>%
  unique %>% 
  cat(sep = "\n")
