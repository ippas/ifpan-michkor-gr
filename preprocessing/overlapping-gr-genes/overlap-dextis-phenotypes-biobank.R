###############################################
# 1. prepare genes from phenotypes from biobank 
###############################################

dir_path <- "data/prs-models-pan-biobank-uk/"

# prepare genes list for phenotypes
phenotypes_biobank_genes_list <- dir_path %>%
  # List all files in the directory with full names
  list.files(full.names = TRUE) %>%
  # Filter files based on the presence of "1e-08" in their names
  keep(~ grepl("1e-08", .x)) %>%
  # Set the names of the list elements to the file names without the extension
  set_names(map(., ~ basename(.x) %>% tools::file_path_sans_ext())) %>%
  # Read YAML content from each file
  map(read_yaml) %>%
  # Extract genes information and name the list elements after the files
  map(~ .x$description$genes)

tissues_clusters <- c("pmid:NA_adrenal-cortex_NA", "pmid:NA_anterior-thigh_NA", "pmid:NA_hypothalamus_NA", 
                      "pmid:NA_kidneys_NA", "pmid:NA_liver_NA", "pmid:NA_lung_NA", 
                      "pmid:NA_perigonadal-adipose-tissue_NA", "pmid:NA_pituitary-gland_NA", "pmid:NA_spleen_NA",
                      "marpiech_tissues_dex_1", "marpiech_tissues_dex_10", "marpiech_tissues_dex_11",
                      "marpiech_tissues_dex_12", "marpiech_tissues_dex_13", "marpiech_tissues_dex_14",
                      "marpiech_tissues_dex_15", "marpiech_tissues_dex_16", "marpiech_tissues_dex_17",
                      "marpiech_tissues_dex_18", "marpiech_tissues_dex_2", "marpiech_tissues_dex_3",
                      "marpiech_tissues_dex_4", "marpiech_tissues_dex_5", "marpiech_tissues_dex_6",
                      "marpiech_tissues_dex_7", "marpiech_tissues_dex_8", "marpiech_tissues_dex_9")

############################################
# 2. calculate chi2 tests for all phenotypes 
############################################
chi2_results_phenotypes <- perform_chi2_tests(c(papers_gene_list[tissues_clusters], phenotypes_biobank_genes_list), hgnc_symbols_vector_v110)

###

###
processing_overlap_results(data = chi2_results_phenotypes,
                           rows_to_filter = !rownames(chi2_results_phenotypes$p_value_matrix) %in% tissues_clusters,
                           cols_to_filter = tissues_clusters[10:27],
                           genes_list = c(papers_gene_list[tissues_clusters], phenotypes_biobank_genes_list)) -> clusters_phenotypes_data


draw_custom_heatmap(
  clusters_phenotypes_data,
  data_type = "significant_uniq_data",
  col_mapping_vector =  clusters_mapping,
  # row_mapping_vector =  setNames(papers_data_preprocessing$label2, papers_data_preprocessing$label),
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
  col_significant = T
)


processing_overlap_results(data = chi2_results_phenotypes,
                           rows_to_filter = !rownames(chi2_results_phenotypes$p_value_matrix) %in% tissues_clusters,
                           cols_to_filter = tissues_clusters[1:9],
                           genes_list = c(papers_gene_list[tissues_clusters], phenotypes_biobank_genes_list)) -> tissue_phenotypes_data


draw_custom_heatmap(
  tissue_phenotypes_data,
  data_type = "significant_uniq_data",
  col_mapping_vector =  tissues_mapping,
  # row_mapping_vector =  setNames(papers_data_preprocessing$label2, papers_data_preprocessing$label),
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
  col_significant = T
)

