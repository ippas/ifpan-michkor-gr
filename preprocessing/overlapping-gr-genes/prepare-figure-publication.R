# Creating a vector in R
papers_mapping <- c("michkor-cells_NA_astrocyte" = "astrocytes, PMID:28381250", 
                 "pmid:18489715_hypothalamus_NA" = "hypothalamus PMID:18489715", 
                 "pmid:19801529_lung_A549" = "lung PMID:19801529",
                 "pmid:21846803_hippocampus_dentate-gyrus" = "dentate gyrus PMID:21846803", 
                 "pmid:22778218_hippocampus_dentate-gyrus" = "dentate gyrus PMID:22778218", 
                 "pmid:23110767_cortex_astrocytes" = "cortex, astrocytes PMID:23110767",
                 "pmid:23339081_striatum_astrocytes" = "striatum, astrocytes PMID:23339081", 
                 "pmid:23339081_striatum_neurones" = "striatum, neurons PMID:23339081", 
                 "pmid:23633533_hippocampus_dentate-gyrus" = "dentate gyrus PMID:23633533",
                 "pmid:23736296_hippocampus_ca1" = "hippocampus, CA1 PMID:23736296", 
                 "pmid:24147833_brain_astrocytes" = "astrocytes PMID:24147833", 
                 "pmid:24147833_brain_mglia" = "microglia PMID:24147833",
                 "pmid:24342991_hippocampus_NA" = "hippocampus PMID:24342991", 
                 "pmid:24777604_embryos_hypothalamic-region_NPSCs" = "embryos, hypothalamic region PMID:24777604", 
                 "pmid:24926665_lung_ASM" = "lung, ASM PMID:24926665",
                 "pmid:26606517_embryos_hypothalamic-region_NPSCs" = "embryos, hypothalamic region PMID:26606517", 
                 "pmid:28500512_prefrontal-cortex_NA" = "prefronal cortex PMID:28500512", 
                 "pmid:31176677_ARC_glia" = "ARC, glia PMID:31176677",
                 "pmid:34362910_hippocampus_NA" = "hippocampus PMID:34362910", 
                 "pmid:35904237_adrenal-gland_Y-1" = "adrenal gland, Y1 PMID:35904237", 
                 "pmid:35948369_embryos_pe" = "embryos, primitive endoderm PMID:35948369")

svg("heatmap_papers.svg", width = 15, height = 11)
draw_custom_heatmap(
  clusters_papers_data,
  data_type = "significant_uniq_data",
  col_mapping_vector =  clusters_mapping,
  row_mapping_vector = papers_mapping,
  # row_mapping_vector =  setNames(papers_data_preprocessing$label2, papers_data_preprocessing$label),
  fdr_threshold = 0.1,
  fdr_thresholds = c(0.01, 0.0001),
  color_rects =  c("green", "#FF00FF"),
  color_rect = "green",
  lwd_rect = 3,
  alpha_rect = 1,
  apply_filling = F,
  color_filling = "gray",
  alpha_filling = 0.6,
  size_filling = 1,
  pch_filling = 16,
  col_significant = T,
  row_label_size = 14,
  col_label_size = 14
)

# Close the SVG device
dev.off()

browseURL("heatmap_papers.svg")


# Open the SVG device
svg("heatmap_output.svg", width = 6, height = 7.2)  # Adjust width and height as needed

  # Draw the heatmap
draw_custom_heatmap(
  clusters_metabolism_data,
  data_type = "significant_uniq_data",
  col_mapping_vector = clusters_mapping,
  row_mapping_vector = metabolism_mapping_vector,
  fdr_threshold = 0.1,
  fdr_thresholds = c(0.05, 0.0001),
  color_rects = c("green", "#FF00FF"),
  color_rect = "green",
  lwd_rect = 3,
  alpha_rect = 1,
  apply_filling = FALSE,
  color_filling = "gray",
  alpha_filling = 0.6,
  size_filling = 1,
  pch_filling = 16,
  row_label_size = 16,
  col_label_size = 16,
  col_significant = TRUE
)

# Close the SVG device
dev.off()

browseURL("heatmap_output.svg")

# Open the SVG device
svg("heatmap_phenotypes.svg", width = 10, height = 12)  # Adjust width and height as needed

draw_custom_heatmap(
  clusters_phenotypes_data,
  data_type = "significant_uniq_data",
  col_mapping_vector =  clusters_mapping,
  row_mapping_vector =  phenotypes_mapping_vector,
  fdr_threshold = 0.1,
  fdr_thresholds = c(0.05, 0.0001),
  color_rects =  c("green", "#FF00FF"),
  color_rect = "green",
  lwd_rect = 3,
  alpha_rect = 1,
  apply_filling = F,
  color_filling = "gray",
  alpha_filling = 0.6,
  size_filling = 1,
  pch_filling = 16,
  col_significant = T,
  row_label_size = 16,
  col_label_size = 16
)

dev.off()

browseURL("heatmap_phenotypes.svg")



##### to suplement
# Open the SVG device
svg("heatmap_phenotypes-0.05.svg", width = 12, height = 14)  
draw_custom_heatmap(
  clusters_phenotypes_data,
  data_type = "significant_uniq_data",
  col_mapping_vector =  clusters_mapping,
  row_mapping_vector =  c(phenotypes_mapping_vector, excluded_vector),
  fdr_threshold = 0.1,
  fdr_thresholds = c(0.05, 0.0001),
  color_rects =  c("green", "#FF00FF"),
  color_rect = "green",
  lwd_rect = 3,
  alpha_rect = 1,
  apply_filling = F,
  color_filling = "gray",
  alpha_filling = 0.6,
  size_filling = 1,
  pch_filling = 16,
  row_label_size = 16,
  col_label_size = 16,
  col_signifi = T
)


dev.off()

browseURL("heatmap_phenotypes-0.05.svg")



###
processing_overlap_results(data = chi2_results_metabolism  ,
                           rows_to_filter = !rownames(chi2_results_metabolism$p_value_matrix) %in% c(tissues_clusters, c("X-11564", "X-11261", "X-21470", "X-21467")),
                           cols_to_filter = tissues_clusters[10:27],
                           genes_list =  metabolism_gene_list,
                           fdr_threshold = 0.05,
                           overlap_threshold = 3) -> clusters_metabolism_data


# tmp$gene_list_sizes <- c(papers_gene_list, genes_list["master_gr_weak"]) %>% sapply(., length)
# Open the SVG device
svg("heatmap_metabolism-0.05.svg", width = 6.5, height = 8)  # Adjust width and height as needed

draw_custom_heatmap(
  clusters_metabolism_data,
  data_type = "significant_uniq_data",
  col_mapping_vector =  clusters_mapping,
  row_mapping_vector = metabolism_mapping_vector,
  fdr_threshold = 0.1,
  fdr_thresholds = c(0.05, 0.0001),
  color_rects =  c("green", "#FF00FF"),
  color_rect = "green",
  lwd_rect = 3,
  alpha_rect = 1,
  apply_filling = F,
  color_filling = "gray",
  alpha_filling = 0.6,
  size_filling = 1,
  pch_filling = 16,
  row_label_size = 16,
  col_label_size = 16,
  col_significant = T,
)

dev.off()

browseURL("heatmap_metabolism-0.05.svg")
