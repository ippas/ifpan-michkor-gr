results_overlap$raw_list$mcluster$significant_uniq_data$rows

results_overlap$raw_list$mcluster$significant_data$df %>% head

results_overlap$raw_list$mcluster$significant_data$df %>% filter(Var2 == "cluster_1")

analyze_gene_list_overlap(
  row_lists = gr_database_blocked_gene_lists$marpiech_cluster_dex_letters,
  col_lists = set_data$gr_list$refine_gene_lists$gene_lists,
  reference_hgnc_vector = hgnc_symbols_vector_v110, 
  keep_original_data = TRUE, 
  fdr_threshold = 0.01, 
  overlap_threshold = 3
) -> results_overlap$raw_list$mclusters_vs_papers_v16 # Można odtworzyć, na zasadzie wczytania każdej wersji bazy i zapisania wyników, ale to zobaczymy co władcy powiedzą, nie filtruj wyników, później
# problem z wynikami: trzeba się zastanowić czy wszystkie listy są potrzebne, raczej dodatkowe ręczne usunięcie
# usunięcię wyników gdzie było dużo punktów czasowych
# usunięcie wyników gdzie nie był badany wpływ GC 
# możliwe, że w pracach gdzie badano różne dawki to można usunąć więcej niż jedną dawkę (ale to już mniej ważne)
# ten etap można wykonać poprzez przeglądnięcie nazw list

results_overlap$raw_list$mclusters_vs_papers_v16

results_overlap$raw_list$mclusters_vs_papers_v16$significant_uniq_data$df %>% 
  .$overlap_genes %>% 
  lapply(., convert_genes_to_vector, split = ",") %>% 
  unlist() %>% 
  table




clusters_papers_data$significant_uniq_data$df %>% filter(Var2 == "cluster_1")

results_overlap$raw_list$mcluster$significant_uniq_data$df %>% filter(Var2 == "cluster_1")

analyze_gene_list_overlap(
  row_lists = gr_database_blocked_gene_lists$marpiech_cluster_dex_letters,
  col_lists = gr_database_blocked_gene_lists$GR_dependent_gene_lists_pmid_tissue_cell_treatment_type_time_dose_regulation$gene_lists,
  reference_hgnc_vector = hgnc_symbols_vector_v110, 
  keep_original_data = TRUE, 
  fdr_threshold = 0.01, 
  overlap_threshold = 3
) -> mclusters_vs_papers_v16



vector_pattern_to_remove <- c("lung-epithelial", "lung-AirwaySmoothMuscle",
                    "brain-organ", "brain-glia",
                    "blood-noMacrophages", "blood-macrophages") %>% 
  paste(., collapse = "|")

master_and_tissue_gene_lists_filtered <- c(gr_database_blocked_gene_lists$master_and_tissue_gene_lists_down, 
                   gr_database_blocked_gene_lists$master_and_tissue_gene_lists_up)

# master_and_tissue_gene_lists_filtered <- master_and_tissue_gene_lists_filtered[!grepl(vector_pattern_to_remove, names(master_and_tissue_gene_lists_filtered))]

master_and_tissue_gene_lists_filtered %>% names

############################
# najpierw wykresy aby wg coś było:
# z porównań wywalić podzielone listy dla brain, blood, lung

# mclusters
analyze_gene_list_overlap(
  row_lists = gr_database_blocked_gene_lists$marpiech_cluster_dex_letters,
  col_lists = master_and_tissue_gene_lists_filtered,
  reference_hgnc_vector = hgnc_symbols_vector_v110, 
  keep_original_data = TRUE, 
  fdr_threshold = 0.01, 
  overlap_threshold = 3
) -> results_overlap_secondary_gene_lists$mapiech_cluster_dex_vs_secondary_gene_lists


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
  overlap_threshold = 3
)


# phenotypes
analyze_gene_list_overlap(
  col_lists = categorized_gene_lists$phenotypes_PanUkBiobank_1e08$gene_lists,
  row_lists = master_and_tissue_gene_lists_filtered,
  reference_hgnc_vector = hgnc_symbols_vector_v110, 
  keep_original_data = TRUE, 
  fdr_threshold = 0.01, 
  overlap_threshold = 3
) -> results_overlap_secondary_gene_lists$phenotypes_1e08_vs_secondary_gene_lists


analyze_gene_list_overlap(
  col_lists = categorized_gene_lists$phenotypes_PanUkBiobank_1e08$gene_lists,
  row_lists =  master_and_tissue_gene_lists_filtered[c("master_up_double_rs","master_down_double_rs")],
  reference_hgnc_vector = hgnc_symbols_vector_v110, 
  keep_original_data = TRUE, 
  fdr_threshold = 0.01, 
  overlap_threshold = 3
) -> results_overlap_secondary_gene_lists$phenotypes_1e08_vs_global_lists

manual_filter_overlap_results_var2 <- function(data, cols_to_remove) {
  # Extracting the 'significant_uniq_data' component from the input data
  sig_uniq_data <- data$manual_filter_overlap_results
  
  print(sig_uniq_data$cols)
  
  # Identifying columns to keep by excluding the columns specified in 'cols_to_remove'
  cols_to_keep <- setdiff(sig_uniq_data$cols, cols_to_remove)
  
  print(cols_to_keep)
  
  # Filtering matrices and dataframe based on 'cols_to_keep'
  filtered_list <- lapply(sig_uniq_data$list, function(matrix) matrix[, cols_to_keep, drop = FALSE])
  filtered_df <- sig_uniq_data$df
  
  # Extracting unique overlap genes from the filtered dataframe
  overlap_genes <- unique(unlist(strsplit(filtered_df$overlap_genes, ",")))
  
  # Creating the manual filter list with updated components
  manual_filter <- list(
    df = filtered_df,
    list = filtered_list,
    rows = sig_uniq_data$rows,  # rows remain the same
    cols = cols_to_keep,
    overlap_genes = overlap_genes
  )
  
  # Adding the manual filter list to the original data
  data$manual_filter_overlap_results <- manual_filter 
  
  # Returning the updated data
  return(data)
}


tmp <- manual_filter_overlap_results(results_overlap_secondary_gene_lists$phenotypes_1e08_vs_secondary_gene_lists, 
                                     rows_to_remove = c("biobankuk-1747-both_sexes-2-hair_colour_natural_before_greying_-EUR",
                                                        "biobankuk-20023-both_sexes--mean_time_to_correctly_identify_matches-EUR",
                                                        "biobankuk-20116-both_sexes-0-smoking_status-EUR",
                                                        "biobankuk-5264-both_sexes--corneal_hysteresis_left_-EUR",
                                                        "biobankuk-5265-both_sexes--corneal_resistance_factor_left_-EUR",
                                                        "biobankuk-5256-both_sexes--corneal_hysteresis_right_-EUR",
                                                        "biobankuk-5264-both_sexes--corneal_hysteresis_left_-EUR",
                                                        "biobankuk-5085-both_sexes--spherical_power_left_-EUR",
                                                        "biobankuk-5255-both_sexes--intra_ocular_pressure_goldmann_correlated_right_-EUR",
                                                        "biobankuk-j45-both_sexes--j45_asthma-EUR",
                                                        "biobankuk-6159-both_sexes-1-pain_type_s_experienced_in_last_month-EUR"))



tmp <- manual_filter_overlap_results_var2(tmp, 
                                     cols_to_remove = c("embryos_up (3)",
                                                        "liver_down (1)",
                                                        "adipose_down (1)",
                                                        "brain-general_down (5)",
                                                        "cartilage_up (1)",
                                                        "kidney_down (3)"))

tmp$manual_filter_overlap_results$df %>% head

# results_overlap_secondary_gene_lists$phenotypes_1e08_vs_secondary_gene_lists$significant_uniq_data$df

svg(width = 17.5, height = 10, "results/google-drive/figures/phenotypes-vs-secondaryList/phenotypes-geneList-filter.svg")

draw_custom_heatmap(
  tmp,
  data_type = "manual_filter_overlap_results",
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

# metabolite traits
analyze_gene_list_overlap(
  col_lists = c(categorized_gene_lists$metabolome_gene_lists$gene_lists, categorized_gene_lists$nightingale_gene_lists$gene_lists),
  row_lists = master_and_tissue_gene_lists_filtered,
  reference_hgnc_vector = hgnc_symbols_vector_v110, 
  keep_original_data = TRUE, 
  fdr_threshold = 0.01, 
  overlap_threshold = 3
) -> results_overlap_secondary_gene_lists$metabolonNightingale_vs_secondary_gene_lists


draw_custom_heatmap(
  results_overlap_secondary_gene_lists$metabolonNightingale_vs_secondary_gene_lists,
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


master_and_tissue_gene_lists_filtered[grepl("brain", names(master_and_tissue_gene_lists_filtered))]

# brain
analyze_gene_list_overlap(
  col_lists =  categorized_gene_lists$phenotypes_PanUkBiobank$gene_lists,
  row_lists = master_and_tissue_gene_lists_filtered[grepl("brain", names(master_and_tissue_gene_lists_filtered))],
  reference_hgnc_vector = hgnc_symbols_vector_v110, 
  keep_original_data = TRUE, 
  fdr_threshold = 0.01, 
  overlap_threshold = 3
) -> results_overlap_secondary_gene_lists$phenotypes_vs_brain_gene_lists


# psychiatric
analyze_gene_list_overlap(
  col_lists =  categorized_gene_lists$phenotypes_PanUkBiobank$gene_lists,
  row_lists = master_and_tissue_gene_lists_filtered,
  reference_hgnc_vector = hgnc_symbols_vector_v110, 
  keep_original_data = TRUE, 
  fdr_threshold = 0.1, 
  overlap_threshold = 1
) -> results_overlap_secondary_gene_lists$phenotypes_vs_secondary_gene_lists

draw_custom_heatmap(
  results_overlap_secondary_gene_lists$phenotypes_vs_brain_gene_lists$significant_uniq_data$df,
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
