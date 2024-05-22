fast_heatmap <- function(data,
                         data_type,
                         overlap_threshold = 3){
  draw_custom_heatmap(
    data,
    data_type = "significant_uniq_data",
    palette = c(
      "pastel_blue"       = "white",
      "pastel_light_blue" = "#f8dedd",
      "white"        = "#f1bcbb",
      "pastel_orange"= "#edacab",
      "pastel_red"   = "#e68a89"
    ),
    # col_mapping_vector =  clusters_mapping,
    # row_mapping_vector = metabolism_mapping_vector,
    # fdr_threshold = 0.01,
    fdr_thresholds = c(0.01, 0.0001),
    color_rects =  c("#4C8D05", "#66023C"),
    color_rect = "green",
    lwd_rect = 2,
    alpha_rect = 1,
    apply_filling = F,
    color_filling = "gray",
    alpha_filling = 0.6,
    size_filling = 1,
    pch_filling = 16,
    col_significant = T,
    row_dend_width = unit(4, "cm"),  # Adjust row dendrogram width
    column_dend_height = unit(3, "cm"),  # Adjust column dendrogram height
    row_names_gp = gpar(fontsize = 16),
    column_names_gp = gpar(fontsize = 16),
    column_names_rot = 45,
    column_names_side = "top",
    overlap_threshold = overlap_threshold
  )
  
}


restore_overlap_analysis <- list()

################################################################################
# mclusters vs papers
analyze_gene_list_overlap(
  row_lists = gr_database_blocked_gene_lists$marpiech_cluster_dex,
  col_lists = gr_database_blocked_gene_lists$GR_dependent_gene_lists_pmid_tissue_cell_treatment_type_regulation$gene_lists,
  reference_hgnc_vector = hgnc_symbols_vector_v110, 
  keep_original_data = TRUE, 
  fdr_threshold = 0.01, 
  overlap_threshold = 3
) -> restore_overlap_analysis$mcluster_vs_papers

gr_database_blocked_gene_lists$GR_dependent_gene_lists_pmid_tissue_cell_treatment_type_regulation$gene_lists %>%
  lapply(., length) %>%
  unname() %>%
  unlist() %>%
  hist(xlab = "Length of Gene Lists", ylab = "Frequency", main = "Histogram of Gene List Lengths")


fast_heatmap(data = restore_overlap_analysis$mcluster_vs_papers,
             data_type = "significant_uniq_data")



################################################################################
analyze_gene_list_overlap(
  row_lists = gr_database_blocked_gene_lists$marpiech_cluster_dex,
  col_lists = categorized_gene_lists$kegg_pathways_gene_lists$gene_lists,
  reference_hgnc_vector = hgnc_symbols_vector_v110, 
  keep_original_data = TRUE, 
  fdr_threshold = 0.01, 
  overlap_threshold = 3
) -> restore_overlap_analysis$papers_vs_kegg_pathways


fast_heatmap(data = restore_overlap_analysis$papers_vs_kegg_pathways,
             data_type = "significant_uniq_data")


################################################################################
secondary_gene_lists_all$rare_gene_lists
