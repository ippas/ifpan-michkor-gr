detailed_tissue_master_df %>% 
  filter(regulation == "up", 
         rank_criterion == "log2ratio",
         size_list == 50) %>% 
  ungroup(hgnc_symbol) %>% 
  mutate(tissue = paste0(tissue, "_", regulation)) %>% 
  mutate(tissue = paste0(tissue, " (", n_lists_filt, ")")) %>% 
  select(c(tissue, hgnc_symbol)) %>% 
  as.data.frame() %>% 
  split(.$tissue, .$hgnc_symbol ) %>% lapply(., function(x) {x$hgnc_symbol}) %>%  
  append(list(master_up_double_rs = regulation_master_double_rank_score_lists$up_log2ratio_50$hgnc_symbol)) -> gr_database_blocked_gene_lists$master_and_tissue_gene_lists_up

detailed_tissue_master_df %>% 
  filter(regulation == "down", 
         rank_criterion == "log2ratio",
         size_list == 50) %>% 
  ungroup(hgnc_symbol) %>%
  mutate(tissue = paste0(tissue, "_", regulation)) %>% 
  mutate(tissue = paste0(tissue, " (", n_lists_filt, ")")) %>% 
  select(c(tissue, hgnc_symbol)) %>% 
  as.data.frame() %>% 
  split(.$tissue, .$hgnc_symbol ) %>% lapply(., function(x) {x$hgnc_symbol}) %>%  
  append(list(master_down_double_rs = regulation_master_double_rank_score_lists$down_log2ratio_50$hgnc_symbol)) -> gr_database_blocked_gene_lists$master_and_tissue_gene_lists_down

results_overlap_secondary_gene_lists <- list()


################################################################################
# marpiech cluster dex

analyze_gene_list_overlap(
  row_lists = gr_database_blocked_gene_lists$marpiech_cluster_dex,
  col_lists = c(gr_database_blocked_gene_lists$master_and_tissue_gene_lists_down, gr_database_blocked_gene_lists$master_and_tissue_gene_lists_up),
  reference_hgnc_vector = hgnc_symbols_vector_v110, 
  keep_original_data = TRUE, 
  fdr_threshold = 0.01, 
  overlap_threshold = 2
) -> results_overlap_secondary_gene_lists$mapiech_cluster_dex_vs_secondary_gene_lists

# fast_heatmap(data = results_overlap_secondary_gene_lists$mapiech_cluster_dex_vs_secondary_gene_lists,
#              data_type = "significant_uniq_data")


svg("results/google-drive/figures/secondary-lists-overlap/secondary-lists-cluster-ovelap-uniq.svg", width = 13, height = 18)

draw_custom_heatmap(
  results_overlap_secondary_gene_lists$mapiech_cluster_dex_vs_secondary_gene_lists,
  data_type = "significant_uniq_data",
  palette = c(
    "pastel_blue"       = "white",
    "pastel_light_blue" = "#f8dedd",
    "white"        = "#f1bcbb",
    "pastel_orange"= "#edacab",
    "pastel_red"   = "#e68a89"
  ),
  fdr_threshold = 0.01,
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
  overlap_threshold = 2
)

dev.off()

write_tsv_xlsx(data =  results_overlap_secondary_gene_lists$mapiech_cluster_dex_vs_secondary_gene_lists$significant_uniq_data$df,
               tsv_file = "results/google-drive/overlap/significant-uniq-marpiech-clusters-vs-secondary-gene-lists.tsv")

write_tsv_xlsx(data =  results_overlap_secondary_gene_lists$mapiech_cluster_dex_vs_secondary_gene_lists$significant_data$df,
               tsv_file = "results/google-drive/overlap/significant-marpiech-clusters-vs-secondary-gene-lists.tsv")

################################################################################
# phenotypes

analyze_gene_list_overlap(
  col_lists = categorized_gene_lists$phenotypes_PanUkBiobank_1e08$gene_lists,
  row_lists = c(gr_database_blocked_gene_lists$master_and_tissue_gene_lists_down, gr_database_blocked_gene_lists$master_and_tissue_gene_lists_up),
  reference_hgnc_vector = hgnc_symbols_vector_v110, 
  keep_original_data = TRUE, 
  fdr_threshold = 0.01, 
  overlap_threshold = 3
) -> results_overlap_secondary_gene_lists$phenotypes_1e08_vs_secondary_gene_lists

# fast_heatmap(data = results_overlap_secondary_gene_lists$phenotypes_1e08_vs_secondary_gene_lists,
#              data_type = "significant_uniq_data")

svg("results/google-drive/figures/secondary-lists-overlap/secondary-lists-phenotypes1e08-ovelap-uniq.svg", 
    width = 24, 
    height = 18)

draw_custom_heatmap(
  results_overlap_secondary_gene_lists$phenotypes_1e08_vs_secondary_gene_lists,
  data_type = "significant_uniq_data",
  palette = c(
    "pastel_blue"       = "white",
    "pastel_light_blue" = "#f8dedd",
    "white"        = "#f1bcbb",
    "pastel_orange"= "#edacab",
    "pastel_red"   = "#e68a89"
  ),
  fdr_threshold = 0.01,
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
  overlap_threshold = 3
)

dev.off()

write_tsv_xlsx(data =  results_overlap_secondary_gene_lists$phenotypes_1e08_vs_secondary_gene_lists$significant_uniq_data$df,
               tsv_file = "results/google-drive/overlap/significant-uniq-phenotypes1e08-vs-secondary-gene-lists.tsv")

write_tsv_xlsx(data =  results_overlap_secondary_gene_lists$phenotypes_1e08_vs_secondary_gene_lists$significant_data$df,
               tsv_file = "results/google-drive/overlap/significant-phenotypes1e08-vs-secondary-gene-lists.tsv")
             
################################################################################
# metabolon

analyze_gene_list_overlap(
  col_lists = categorized_gene_lists$metabolome_gene_lists$gene_lists,
  row_lists = c(gr_database_blocked_gene_lists$master_and_tissue_gene_lists_down, gr_database_blocked_gene_lists$master_and_tissue_gene_lists_up),
  reference_hgnc_vector = hgnc_symbols_vector_v110, 
  keep_original_data = TRUE, 
  fdr_threshold = 0.01, 
  overlap_threshold = 2
) -> results_overlap_secondary_gene_lists$metabolon_vs_secondary_gene_lists


# fast_heatmap(data = results_overlap_secondary_gene_lists$metabolon_vs_secondary_gene_lists,
#              data_type = "significant_uniq_data")

svg("results/google-drive/figures/secondary-lists-overlap/secondary-lists-metabolon-ovelap-uniq.svg", 
    width = 14, 
    height = 11)

draw_custom_heatmap(
  results_overlap_secondary_gene_lists$metabolon_vs_secondary_gene_lists,
  data_type = "significant_uniq_data",
  palette = c(
    "pastel_blue"       = "white",
    "pastel_light_blue" = "#f8dedd",
    "white"        = "#f1bcbb",
    "pastel_orange"= "#edacab",
    "pastel_red"   = "#e68a89"
  ),
  fdr_threshold = 0.01,
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
  overlap_threshold = 2
)

dev.off()

write_tsv_xlsx(data =  results_overlap_secondary_gene_lists$metabolon_vs_secondary_gene_lists$significant_uniq_data$df,
               tsv_file = "results/google-drive/overlap/significant-uniq-metabolome-vs-secondary-gene-lists.tsv")

write_tsv_xlsx(data =  results_overlap_secondary_gene_lists$metabolon_vs_secondary_gene_lists$significant_data$df,
               tsv_file = "results/google-drive/overlap/significant-metabolome-vs-secondary-gene-lists.tsv")

################################################################################
# nightingale

analyze_gene_list_overlap(
  col_lists = categorized_gene_lists$nightingale_gene_lists$gene_lists,
  row_lists = c(gr_database_blocked_gene_lists$master_and_tissue_gene_lists_down, gr_database_blocked_gene_lists$master_and_tissue_gene_lists_up),
  reference_hgnc_vector = hgnc_symbols_vector_v110, 
  keep_original_data = TRUE, 
  fdr_threshold = 0.01, 
  overlap_threshold = 2
) -> results_overlap_secondary_gene_lists$nightingale_vs_secondary_gene_lists

# fast_heatmap(data = results_overlap_secondary_gene_lists$nightingale_vs_secondary_gene_lists,
#              data_type = "significant_uniq_data")

svg("results/google-drive/figures/secondary-lists-overlap/secondary-lists-nightingale-ovelap-uniq.svg", 
    width = 11, 
    height = 9)

draw_custom_heatmap(
  results_overlap_secondary_gene_lists$nightingale_vs_secondary_gene_lists,
  data_type = "significant_uniq_data",
  palette = c(
    "pastel_blue"       = "white",
    "pastel_light_blue" = "#f8dedd",
    "white"        = "#f1bcbb",
    "pastel_orange"= "#edacab",
    "pastel_red"   = "#e68a89"
  ),
  fdr_threshold = 0.01,
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
  overlap_threshold = 2
)

dev.off()

write_tsv_xlsx(data =  results_overlap_secondary_gene_lists$nightingale_vs_secondary_gene_lists$significant_uniq_data$df,
               tsv_file = "results/google-drive/overlap/significant-uniq-nightingale-vs-secondary-gene-lists.tsv")

write_tsv_xlsx(data =  results_overlap_secondary_gene_lists$nightingale_vs_secondary_gene_lists$significant_data$df,
               tsv_file = "results/google-drive/overlap/significant-nightingale-vs-secondary-gene-lists.tsv")
