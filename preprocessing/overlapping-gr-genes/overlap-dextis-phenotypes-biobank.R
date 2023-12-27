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
# Creating a vector from the provided list
phenotypes_mapping_vector <- c(
  "biobankuk-101270-both_sexes--other_bread_intake-EUR-1e-08" = "other bread intake ID:30620",
  "biobankuk-20004-both_sexes-1479-operation_code-EUR-1e-08" = "varicose vein surgery ID:454",
  "biobankuk-30280-both_sexes--immature_reticulocyte_fraction-EUR-1e-08" = "immature reticulocyte fraction ID:30280",
  "biobankuk-30290-both_sexes--high_light_scatter_reticulocyte_percentage-EUR-1e-08" = "reticulocyte percentage ID:30290",
  "biobankuk-30300-both_sexes--high_light_scatter_reticulocyte_count-EUR-1e-08" = "reticulocyte count ID:30300",
  "biobankuk-30620-both_sexes--alanine_aminotransferase-EUR-1e-08" = "alanine aminotransferase ID:30620",
  "biobankuk-3064-both_sexes--peak_expiratory_flow_pef_-EUR-1e-08" = "peak expiratory flow ID:3064",
  "biobankuk-30870-both_sexes--triglycerides-EUR-1e-08" = "triglycerides ID:30870",
  "biobankuk-4080-both_sexes--systolic_blood_pressure_automated_reading-EUR-1e-08" = "systolic blood pressure ID:4080",
  "biobankuk-454_1-both_sexes--varicose_veins_of_lower_extremity-EUR-1e-08" = "varicose veins of lower extremity ID:131402",
  "biobankuk-454-both_sexes--varicose_veins-EUR-1e-08" = "varicose veins ID:131408",
  "biobankuk-5098-both_sexes--6mm_weak_meridian_right_-EUR-1e-08" = "6mm weak meridian right ID:5098",
  "biobankuk-sbp-both_sexes--systolic_blood_pressure_automated_reading_adjusted_by_medication-EUR-1e-08" = "systolic blood presure adjusted by medication",
  "biobankuk-30240-both_sexes--reticulocyte_percentage-EUR-1e-08" = "reticulocyte percentage ID:30240"
)

# Exclude elements from data_vector that are in phenotypes_mapping_vector
excluded_vector <- setdiff(clusters_phenotypes_data$significant_uniq_data$rows, names(phenotypes_mapping_vector))

excluded_vector <- c("biobankuk-1528-both_sexes--water_intake-EUR-1e-08" = "watere intake ID:1528",
                     "biobankuk-1717-both_sexes--skin_colour-EUR-1e-08" = "skin colour ID:1717",
                     "biobankuk-1727-both_sexes--ease_of_skin_tanning-EUR-1e-08" = "ease of skin tanning ID:1727",
                     "biobankuk-1747-both_sexes-1-hair_colour_natural_before_greying_-EUR-1e-08" = "blond hair colour natural ID:1747",
                     "biobankuk-23118-both_sexes--leg_predicted_mass_left_-EUR-1e-08" = "leg predicted mass left ID:23118",
                     "biobankuk-30120-both_sexes--lymphocyte_count-EUR-1e-08" = "lymphocyte count ID:30120",
                     "biobankuk-30240-both_sexes--reticulocyte_percentage-EUR-1e-08" = "reticulocyte percentage ID:30240",
                     "biobankuk-30630-both_sexes--apolipoprotein_a-EUR-1e-08" = "apolipoprotein A ID:30630",
                     "biobankuk-30660-both_sexes--direct_bilirubin-EUR-1e-08" = "direct bilirubin ID:30660",
                     "biobankuk-30830-both_sexes--shbg-EUR-1e-08" = "shbg ID:30830",
                     "biobankuk-5132-both_sexes--3mm_strong_meridian_right_-EUR-1e-08" = "3mm strong meridian right ID:5132",
                     "biobankuk-5135-both_sexes--3mm_strong_meridian_left_-EUR-1e-08" = "3mm strong meridian left ID:5135",
                     "biobankuk-5256-both_sexes--corneal_hysteresis_right_-EUR-1e-08" = "corneal hysteresis right ID:5256")

c(papers_gene_list[tissues_clusters], phenotypes_biobank_genes_list) %>% 
    lapply(., unique) -> phenotype_list

chi2_results_phenotypes <- perform_chi2_tests(phenotype_list, hgnc_symbols_vector_v110)

###

###
processing_overlap_results(data = chi2_results_phenotypes,
                           rows_to_filter = !rownames(chi2_results_phenotypes$p_value_matrix) %in% tissues_clusters,
                           cols_to_filter = tissues_clusters[10:27],
                           genes_list = c(papers_gene_list[tissues_clusters], phenotypes_biobank_genes_list),
                           fdr_threshold = 0.05,
                           overlap_threshold = 3) -> clusters_phenotypes_data


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
  col_signifi = T
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

