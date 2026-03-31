# ##############################################################################
# ---- prepare data ----
# ##############################################################################

factors_rsidGenes_P1e2locus100kb <- read.table(file = "data/factors-pmid38965376/annotated-rsid-p1e2/annotated_rsid_locus100kbp_p1e2.tsv",
                                               sep = "\t",
                                               header = TRUE) %>% 
  group_by(factor_id, gene_symbol) %>% 
  slice_min(order_by = pvalue, with_ties = FALSE) %>% 
  ungroup


AllBrain2BrainSignatures2Global %>% 
  lapply()

# ##############################################################################
# ---- function ----
# ##############################################################################
process_overlap_heatmap_plot_only <- function(
    factors_df,
    phenotype_col = "factor_name",   # << teraz jasno: kolumna z fenotypem
    type = c("locus", "geneCenter"),
    window_kb = NULL,
    pvalue_threshold = 1e-5,
    drawing_overlap_threshold = 3,
    pvalue_overlap_threshold = 0.05,
    plot_title_prefix = "",
    col_only_n_genes = TRUE,
    color_scale_range = c(0, 4),
    reference_gene_lists = NULL,
    return_df = FALSE
) {
  type <- match.arg(type)
  
  # sprawdzamy czy podana kolumna istnieje
  if (!phenotype_col %in% colnames(factors_df)) {
    stop(sprintf("Column '%s' not found in factors_df", phenotype_col))
  }
  
  # jeśli nie podano własnych list genów, bierzemy lite_grSignatures jako domyślne
  if (is.null(reference_gene_lists)) {
    reference_gene_lists <- split(lite_grSignatures$hgnc_symbol,
                                  lite_grSignatures$signature_name)
  }
  
  # grupowanie po kolumnie fenotypu
  factor_gene_list <- factors_df %>%
    dplyr::filter(pvalue < pvalue_threshold) %>%
    dplyr::group_by(.data[[phenotype_col]]) %>%
    dplyr::summarise(genes = list(unique(gene_symbol)), .groups = "drop") %>%
    tibble::deframe()
  
  gene_lists_all <- c(reference_gene_lists, factor_gene_list)
  
  chi2_results <- perform_chi2_tests(
    gene_lists_all,
    total_genes = hgnc_symbols_vector_v110
  )
  
  rows_to_filter <- names(reference_gene_lists)
  cols_to_filter <- factors_df %>%
    dplyr::filter(pvalue < pvalue_threshold) %>%
    dplyr::pull(.data[[phenotype_col]]) %>%
    unique()
  
  overlap <- processing_overlap_results(
    data = chi2_results,
    rows_to_filter = rows_to_filter,
    cols_to_filter = cols_to_filter,
    overlap_threshold = 0,
    fdr_threshold = 1,
    genes_list = gene_lists_all
  )
  
  cat("\n=== original_data$df (p_value < ", pvalue_overlap_threshold, ") ===\n", sep = "")
  print(overlap$original_data$df %>% dplyr::filter(p_value < pvalue_overlap_threshold))
  cat("\n")
  
  title_str <- sprintf("%s     %s +/- %s",
                       plot_title_prefix, type,
                       ifelse(!is.null(window_kb), paste0(window_kb, "kb"), "?"))
  
  p <- draw_custom_heatmap_ggplot_v3(
    data_list = overlap,
    title = title_str,
    title_size = 22,
    data_type = "original_data",
    p_thresholds = c(0.05, 0.01),
    color_rects = c("#4C8D05", "#66023C"),
    overlap_threshold = drawing_overlap_threshold,
    apply_filling = FALSE,
    col_only_n_genes = col_only_n_genes,
    col_significant = TRUE,
    color_scale_range = color_scale_range,
    palette = c("white", "#f8dedd", "#f1bcbb", "#edacab", "#e68a89"),
    text_color = "black",
    axis_text_angle = 90,
    axis_text_hjust = 0,
    text_size_tile = 4,
    text_size_axis = 16,
    text_size_legend = 18
  )
  
  result <- list(plot = p, overlap = overlap)
  
  if (return_df) {
    result$df <- overlap$original_data$df
  }
  
  return(result)
}

# ##############################################################################
# ---- analysis ----
# ##############################################################################
metaphenotypes_AllBrain2Brain2Global <- process_overlap_heatmap_plot_only(factors_rsidGenes_P1e2locus50kb,
                                        type = "locus",
                                        window_kb = 50,
                                        reference_gene_lists = AllBrain2BrainSignatures2Global,
                                        plot_title_prefix = "A",
                                        pvalue_threshold = 0.00001,
                                        color_scale_range = c(0, 8),
                                        col_only_n_genes = F)

p1$overlap$original_data$df %>% filter(Var2 == "clinical_anxiety_and_depression") %>% 
  filter(overlap_genes != "") %>% 
  .$overlap_genes %>% 
  strsplit(split = ",") %>% 
  unlist %>% unique() %>% 
  cat(sep = "\n")
