# ##############################################################################
# ---- uses data ----
# ##############################################################################
AllBrainGeneDf
AllBrainGeneLists


# ##############################################################################
# ---- functions ----
# ##############################################################################
plot_full_overlap_heatmap <- function(
    gene_lists,
    total_genes = hgnc_symbols_vector_v110,
    plot_title = "Full overlap heatmap",
    color_scale_range = c(0, 4),
    drawing_overlap_threshold = 3,
    col_only_n_genes = TRUE,
    return_df = FALSE,
    verbose = TRUE
) {
  
  #' Plot full overlap heatmap (all-vs-all)
  #'
  #' This function computes and visualizes gene overlaps between *all* provided
  #' gene lists using chi-square statistics, without any filtering.  
  #' The result is a **square heatmap**, where both rows and columns represent
  #' the same set of gene lists.
  #'
  #' @param gene_lists Named list of character vectors, each containing gene symbols.
  #' @param total_genes Character vector of background gene symbols (default: hgnc_symbols_vector_v110).
  #' @param plot_title Title for the heatmap.
  #' @param color_scale_range Numeric vector of length 2 defining the scale range for color intensity.
  #' @param drawing_overlap_threshold Minimal overlap (n genes) to be shown as non-empty.
  #' @param col_only_n_genes Logical; if TRUE, tiles are colored only by number of overlapping genes.
  #' @param return_df Logical; if TRUE, also returns data frame with overlap results.
  #' @param verbose Logical; print progress messages.
  #'
  #' @return A list with:
  #' \item{plot}{ggplot2 object — the generated heatmap.}
  #' \item{overlap}{List with overlap matrices and stats.}
  #' \item{df}{(optional) Data frame with overlap statistics, if return_df = TRUE.}
  #'
  #' @examples
  #' plot_full_overlap_heatmap(list(
  #'   ListA = c("FKBP5", "NR3C1", "PER1"),
  #'   ListB = c("FKBP5", "TSC22D3", "SGK1"),
  #'   ListC = c("PER1", "SGK1", "FKBP4")
  #' ))
  #'
  if (!is.list(gene_lists)) stop("❌ gene_lists must be a named list of gene vectors")
  if (is.null(names(gene_lists))) stop("❌ gene_lists must have names for labeling the heatmap")
  
  if (verbose) message("▶ Number of gene lists provided: ", length(gene_lists))
  
  # 1️⃣ Perform chi2 tests on all combinations
  chi2_results <- perform_chi2_tests(
    gene_lists,
    total_genes = total_genes
  )
  
  # 2️⃣ Ensure consistency and unique naming
  rown <- rownames(chi2_results$p_value_matrix)
  coln <- colnames(chi2_results$p_value_matrix)
  rownames(chi2_results$p_value_matrix) <- make.unique(rown)
  colnames(chi2_results$p_value_matrix) <- make.unique(coln)
  rownames(chi2_results$number_overlap_matrix) <- make.unique(rown)
  colnames(chi2_results$number_overlap_matrix) <- make.unique(coln)
  rownames(chi2_results$chi2_value_matrix) <- make.unique(rown)
  colnames(chi2_results$chi2_value_matrix) <- make.unique(coln)
  rownames(chi2_results$overlap_genes_matrix) <- make.unique(rown)
  colnames(chi2_results$overlap_genes_matrix) <- make.unique(coln)
  
  # 3️⃣ Compute overlap results (no filtering, all-vs-all)
  overlap <- processing_overlap_results(
    data = chi2_results,
    rows_to_filter = rownames(chi2_results$p_value_matrix),
    cols_to_filter = colnames(chi2_results$p_value_matrix),
    overlap_threshold = 0,
    fdr_threshold = 1,
    genes_list = gene_lists
  )
  
  if (verbose)
    message("▶ Overlap results: ", nrow(overlap$original_data$df), " total comparisons")
  
  # 4️⃣ Draw square heatmap
  p <- draw_custom_heatmap_ggplot_v3(
    data_list = overlap,
    title = plot_title,
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
  if (return_df) result$df <- overlap$original_data$df
  return(result)
}



# Wywołanie funkcji
tmp_ <- plot_full_overlap_heatmap(
  gene_lists = AllBrainGeneLists,
  total_genes = hgnc_symbols_vector_v110,   # referencyjny zbiór genów (możesz podać własny)
  plot_title = "Full overlap of GR-dependent gene sets",
  color_scale_range = c(0, 4),
  drawing_overlap_threshold = 2,
  col_only_n_genes = TRUE,
  return_df = TRUE,
  verbose = TRUE
)

tmp$overlap$original_data$df %>% select(-fdr) %>% 
  group_by(Var1) %>% 
  nest() %>% 
  mutate(data = map(data, ~ .x %>% mutate(fdr = p.adjust(p_value, method = "fdr")))) %>% 
  mutate(n_signif_association = map(data, ~ .x %>% filter(fdr < 0.1) %>% nrow())) %>% 
  unnest(n_signif_association) %>% 
  filter(n_signif_association < 10) %>% .$Var1


