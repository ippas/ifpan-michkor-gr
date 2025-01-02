#############
# functions #
#############
create_mapping_vector <- function(input_vector) {
  # Extract the last element after the underscore
  extracted_elements <- sapply(strsplit(input_vector, "_"), tail, 1)
  
  # Format the extracted elements
  # formatted_elements <- paste("cluster", extracted_elements)
  
  # Create the named vector
  mapping_vector <- setNames(extracted_elements, input_vector)
  
  return(mapping_vector)
}

#######################
# metabolism analysis #
#######################
# prepare data to chi2
gr_gene_database_preproccesing %>% 
  filter(
    grepl("omicspred_metabolon", source) |
      grepl("omicspred_nithingale_", source) |
      grepl("marpiech_tissues_dex", source) |
      gene_list_index == "all_significant_genes_marpiech"
  ) %>% 
  mutate(source = ifelse(source == "pmid:NA", "marpiech_tissues_dex", source)) %>%
  extract_keys_values("info", c("BiomarkerName", "BiochemicalName")) %>%
  mutate(label = ifelse(
    is.na(BiochemicalName),
    ifelse(is.na(BiomarkerName), gene_list_index, BiomarkerName),
    BiochemicalName
  )) %>%
  mutate(label = ifelse(
    source == "marpiech_tissues_dex",
    paste0(source, "_", gene_list_number),
    label
  )) %>% 
  mutate(label = ifelse(
    gene_list_index == "all_significant_genes_marpiech",
    tissue,
    label
  ))  %>%
  select(c(source, gene_list_index, gene_list_number, hgnc_symbol, label)) %>%
  filter(!is.na(hgnc_symbol)) %>% 
  distinct() %>%
  filter(!(hgnc_symbol %in% hgnc_to_remove)) %>% 
  group_by(gene_list_number) %>%
  nest() %>%
  mutate(n_genes = map_int(data, ~ nrow(.))) %>% 
  unnest(cols = c(data)) %>% 
  ungroup() %>% 
  filter(n_genes > 2) -> metabolism_data_preprocessing

split(
  metabolism_data_preprocessing$hgnc_symbol,
  metabolism_data_preprocessing$label
) -> metabolism_gene_list

lapply(metabolism_gene_list, unique) -> metabolism_gene_list

# calculate chi2 tests
chi2_results_metabolism 


# prepare mapping vector to figure
tissues_clusters <- c("adrenal-cortex", "anterior-thigh", "hypothalamus", 
                      "kidneys", "liver", "lung", 
                      "perigonadal-adipose-tissue", "pituitary-gland", "spleen",
                      "marpiech_tissues_dex_1", "marpiech_tissues_dex_10", "marpiech_tissues_dex_11",
                      "marpiech_tissues_dex_12", "marpiech_tissues_dex_13", "marpiech_tissues_dex_14",
                      "marpiech_tissues_dex_15", "marpiech_tissues_dex_16", "marpiech_tissues_dex_17",
                      "marpiech_tissues_dex_18", "marpiech_tissues_dex_2", "marpiech_tissues_dex_3",
                      "marpiech_tissues_dex_4", "marpiech_tissues_dex_5", "marpiech_tissues_dex_6",
                      "marpiech_tissues_dex_7", "marpiech_tissues_dex_8", "marpiech_tissues_dex_9")

metabolism_mapping_vector <- c("1-oleoyl-3-linoleoyl-glycerol(18" = "1-oleoyl-3-linoleoyl-glycerol",
                               "3-methylglutarylcarnitine(2)" = "3-methylglutarylcarnitine",
                               "4-androsten-3beta17beta-diolmonosulfate(2)" = "4-androsten-3beta 17beta-diol monosulfate",
                               "TriglyceridesinHDL" = "HDL triglycerides",
                               "AveragediameterforHDLparticles" = "HDL diameter",
                               "AveragediameterforLDLparticles" = "LDL diameter",
                               "betaine" = "betaine",
                               "HDLcholesterol" = "HDL cholesterol",
                               "CholesterolinmediumHDL" = "cholesterol medium HDL",
                               "CholesterolinsmallHDL" = "Small HDL cholesterol",
                               "TotallipidsinsmallHDL" = "total lipids in small HDL",
                               "Degreeofunsaturation" = "degree of unsaturation",
                               "Glutamine" = "glutamine",
                               "TriglyceridesinverylargeHDL" = "Very large HDL triglycerides ",
                               "stearoylsphingomyelin(d18" = "stearoyl Sphingomyelin",
                               "TriglyceridesinLDL" = "LDL triglycerides",
                               "TriglyceridesinmediumHDL" = "Medium HDL triglycerides",
                               "Tyrosine" = "tyrosine")


processing_overlap_results(data = chi2_results_metabolism  ,
                           rows_to_filter = !rownames(chi2_results_metabolism$p_value_matrix) %in% c(tissues_clusters),
                           cols_to_filter = tissues_clusters[10:27],
                           genes_list =  metabolism_gene_list,
                           fdr_threshold = 0.01,
                           overlap_threshold = 3) -> clusters_metabolism_data


draw_custom_heatmap(
  clusters_metabolism_data,
  data_type = "significant_uniq_data",
  col_mapping_vector =  clusters_mapping,
  row_mapping_vector = metabolism_mapping_vector,
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
)


processing_overlap_results(data = chi2_results_metabolism ,
                           rows_to_filter = !rownames(chi2_results_metabolism$p_value_matrix) %in% tissues_clusters,
                           cols_to_filter = tissues_clusters[1:9],
                           genes_list = metabolism_gene_list) -> tissue_metabolism_data

draw_custom_heatmap(
  tissue_metabolism_data,
  data_type = "significant_data",
  fdr_threshold = 0.1,
  fdr_thresholds = c(0.05, 0.0001),
  color_rects =  c("green", "#FF00FF"),
  color_rect = "green",
  lwd_rect = 3,
  alpha_rect = 1,
  apply_filling = F,
  color_filling = "green",
  alpha_filling = 0.6,
  size_filling = 1,
  pch_filling = 16,
  col_significant = T
) -> p1
