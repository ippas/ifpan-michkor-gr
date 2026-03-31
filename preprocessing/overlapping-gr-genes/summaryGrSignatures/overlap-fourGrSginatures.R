
lite_grSignatures %>% 
  filter(grepl("brain",signature_name)) %>% .$hgnc_symbol %>% unique

gene_lists_all <- split(lite_grSignatures$hgnc_symbol, lite_grSignatures$signature_name)


chi2_results <- perform_chi2_tests(gene_lists_all, total_genes = hgnc_symbols_vector_v110)

rows_to_filter <- unique(lite_grSignatures$signature_name)

overlap <- processing_overlap_results(
  data = chi2_results,
  rows_to_filter = rows_to_filter,
  cols_to_filter = rows_to_filter,
  overlap_threshold = 0,
  fdr_threshold = 1,
  genes_list = gene_lists_all
)

cat("\n=== original_data$df (p_value < ", pvalue_overlap_threshold, ") ===\n", sep = "")
print(overlap$original_data$df %>% dplyr::filter(p_value < pvalue_overlap_threshold))
cat("\n")

title_str <- sprintf("%s     %s +/- %s", plot_title_prefix, type,
                     ifelse(!is.null(window_kb), paste0(window_kb, "kb"), "?"))

dev.off()

svg("data/summary-signatures/figures/overlap-grFourSignatures.svg", width = 7, height = 7)
draw_custom_heatmap_ggplot_v3(
  data_list = overlap,
  # title = title_str,
  title_size = 22,
  data_type = "original_data",
  p_thresholds = c(0.05, 0.01),
  color_rects = c("#4C8D05", "#66023C"),
  overlap_threshold = drawing_overlap_threshold,
  apply_filling = FALSE,
  # col_only_n_genes = col_only_n_genes,
  col_significant = F,
  color_scale_range = c(0,12),
  palette = c("white", "#f8dedd", "#f1bcbb", "#edacab", "#e68a89"),
  text_color = "black",
  axis_text_angle = 45,
  axis_text_hjust = 0,
  text_size_tile = 4,
  text_size_axis = 16,
  text_size_legend = 18
)

dev.off()
