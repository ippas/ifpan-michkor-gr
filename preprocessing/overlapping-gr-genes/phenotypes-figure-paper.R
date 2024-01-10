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
# 
# ############################################
# # 2. calculate chi2 tests for all phenotypes 
# ############################################
c(papers_gene_list[tissues_clusters[10:27]], phenotypes_biobank_genes_list) %>% 
  lapply(., unique) -> phenotypes_clusters_list 


chi2_results_phenotypes <- perform_chi2_tests(phenotypes_clusters_list, hgnc_symbols_vector_v110)

###

###
processing_overlap_results(data = chi2_results_phenotypes,
                           rows_to_filter = !rownames(chi2_results_phenotypes$p_value_matrix) %in% tissues_clusters,
                           # rows_to_filter = tissues_clusters[10:27],
                           cols_to_filter = tissues_clusters[10:27],
                           overlap_threshold = 3,
                           fdr_threshold = 0.01,
                           genes_list = phenotypes_clusters_list) -> clusters_phenotypes_data



clusters_phenotypes_data$significant_uniq_data$df %>% 
  filter(fdr > 0.01) %>% .$Var1 %>% unique() -> phenotypes_to_remove


manual_filter_overlap_results(clusters_phenotypes_data, rows_to_remove = phenotypes_to_remove) -> clusters_papers_data


phenotyepes_mapping <- c(
    "biobankuk-101270-both_sexes--other_bread_intake-EUR-1e-08" = "other bread intake ID:30620",
    "biobankuk-20004-both_sexes-1479-operation_code-EUR-1e-08" = "varicose vein surgery ID:1479",
    "biobankuk-30240-both_sexes--reticulocyte_percentage-EUR-1e-08" = "reticulocyte percentage ID:30240",
    "biobankuk-30280-both_sexes--immature_reticulocyte_fraction-EUR-1e-08" = "immature reticulocyte fraction ID:30280",
    "biobankuk-30290-both_sexes--high_light_scatter_reticulocyte_percentage-EUR-1e-08" = "reticulocyte percentage ID:30290" ,
    "biobankuk-30300-both_sexes--high_light_scatter_reticulocyte_count-EUR-1e-08" = "reticulocyte count ID:30300",
    "biobankuk-30620-both_sexes--alanine_aminotransferase-EUR-1e-08" = "alanine aminotransferase ID:30620",
    "biobankuk-3064-both_sexes--peak_expiratory_flow_pef_-EUR-1e-08" = "peak expiratory flow ID:3064",
    "biobankuk-30870-both_sexes--triglycerides-EUR-1e-08" = "triglycerides ID:30870",
    "biobankuk-4080-both_sexes--systolic_blood_pressure_automated_reading-EUR-1e-08" = "systolic blood pressure ID:4080",
    "biobankuk-454_1-both_sexes--varicose_veins_of_lower_extremity-EUR-1e-08" = "varicose veins of lower extremity ID:131402",
    "biobankuk-454-both_sexes--varicose_veins-EUR-1e-08" = "varicose veins ID:131408",
    "biobankuk-5098-both_sexes--6mm_weak_meridian_right_-EUR-1e-08" = "6mm weak meridian right ID:5098",
    "biobankuk-sbp-both_sexes--systolic_blood_pressure_automated_reading_adjusted_by_medication-EUR-1e-08" = "systolic blood presure adjusted by medication"
  )

# Replacing " ID" with "                               ID"
phenotypes_mapping <- gsub(" ID", "                               ID", phenotyepes_mapping)


clusters_papers_data$significant_uniq_data$df


# Printing the modified vector
print(phenotyepes_mapping)

draw_custom_heatmap(
  clusters_phenotypes_data,
  data_type = "significant_uniq_data",
  palette = c(
    "pastel_blue"       = "white",
    "pastel_light_blue" = "#f8dedd",
    "white"        = "#f1bcbb",
    "pastel_orange"= "#edacab",
    "pastel_red"   = "#e68a89"
  ),
  col_mapping_vector =  clusters_mapping,
  row_mapping_vector =  phenotyepes_mapping,
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
svg("results/figures/figure-phenotype-metabolome/phenotypes-clusters.svg", width = 13, height = 9)

draw_custom_heatmap(
  clusters_phenotypes_data,
  data_type = "significant_uniq_data",
  palette = c(
    "pastel_blue"       = "white",
    "pastel_light_blue" = "#f8dedd",
    "white"        = "#f1bcbb",
    "pastel_orange"= "#edacab",
    "pastel_red"   = "#e68a89"
  ),
  col_mapping_vector =  clusters_mapping,
  row_mapping_vector =  phenotyepes_mapping,
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

