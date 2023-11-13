
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

# calculate chi2 tests
chi2_results_metabolism <- perform_chi2_tests(metabolism_gene_list, hgnc_symbols_vector_v110)

tissues_clusters <- c("adrenal-cortex", "anterior-thigh", "hypothalamus", 
                      "kidneys", "liver", "lung", 
                      "perigonadal-adipose-tissue", "pituitary-gland", "spleen",
                      "marpiech_tissues_dex_1", "marpiech_tissues_dex_10", "marpiech_tissues_dex_11",
                      "marpiech_tissues_dex_12", "marpiech_tissues_dex_13", "marpiech_tissues_dex_14",
                      "marpiech_tissues_dex_15", "marpiech_tissues_dex_16", "marpiech_tissues_dex_17",
                      "marpiech_tissues_dex_18", "marpiech_tissues_dex_2", "marpiech_tissues_dex_3",
                      "marpiech_tissues_dex_4", "marpiech_tissues_dex_5", "marpiech_tissues_dex_6",
                      "marpiech_tissues_dex_7", "marpiech_tissues_dex_8", "marpiech_tissues_dex_9")


################################################################################
###########################################
# 1. clusters vs metabolism visualization #
###########################################
# Extract significant results for clusters vs metabolism
# This block extracts data from chi-squared results related to metabolism, excluding rows that match tissue clusters.
# It then renames the variables for clarity, adjusts p-values for multiple testing using the False Discovery Rate (FDR) method,
# sorts the results by FDR, and filters for significant results (FDR < 0.05) with more than one overlapping element.
extract_data(chi2_results_metabolism,  
             !rownames(chi2_results_metabolism$p_value_matrix) %in% tissues_clusters,
             tissues_clusters[10:27]
) %>% 
  rename(cluster = Var1, metabolism = Var2) %>% 
  mutate(fdr = p.adjust(p_value, method = "fdr")) %>% 
  arrange(fdr) %>% 
  filter(fdr < 0.05) %>% 
  filter(number_overlap > 1)

# Prepare genes overlapping with metabolism
# This block repeats the extraction process and focuses on overlapping genes.
# It splits the genes by comma, flattens the list, ensures gene uniqueness, and prints each gene on a new line.
extract_data(chi2_results_metabolism,  
             !rownames(chi2_results_metabolism$p_value_matrix) %in% tissues_clusters,
             tissues_clusters[10:27]) %>%
  rename(metabolism = Var1, cluster = Var2) %>% 
  mutate(fdr = p.adjust(p_value, method = "fdr")) %>% 
  arrange(fdr) %>% 
  filter(fdr < 0.05) %>% 
  filter(number_overlap > 1) %>% 
  .$overlap_genes %>% 
  strsplit(., split = ",") %>% 
  unlist %>% 
  unique() %>% 
  walk(., ~cat(.x, "\n"))

# Create a named vector for cluster to metabolism mapping
# This vector is used for labeling in the heatmap later.
cluster_mapping_metabolism <- c("marpiech_tissues_dex_9" = "cluster 9", 
                                "marpiech_tissues_dex_12" = "cluster 12", 
                                "marpiech_tissues_dex_18" = "upregulated genes")

# Visualization of unique significant gene list
# This block prepares a filtered dataset to identify unique significant genes.
# It groups by overlapping genes and clusters, selects the top entry after arranging by FDR,
# and creates a mapping vector for metabolism to clusters.
extract_data(chi2_results_metabolism,  
             tissues_clusters[10:27],
             !rownames(chi2_results_metabolism$p_value_matrix) %in% tissues_clusters) %>%
  rename(cluster = Var1, metabolism = Var2) %>% 
  mutate(fdr = p.adjust(p_value, method = "fdr")) %>% 
  arrange(fdr) %>% 
  filter(fdr < 0.05) %>% 
  filter(number_overlap > 1) %>% 
  filter(!(metabolism %in% c("X-11564", "X-11261", "X-21470", "X-21467"))) %>%
  group_by(overlap_genes, cluster) %>% 
  nest() %>% 
  mutate(data = map(data, ~ .x %>% 
                      arrange(fdr) %>% 
                      slice(1))) %>% 
  unnest() %>% 
  ungroup %>% 
  as.data.frame() %>% 
  .$metabolism %>% 
  as.character() %>%
  unique() %>% 
  create_mapping_vector() -> metabolism_clusters_mapping

# Generate a heatmap for metabolism clusters
# This heatmap visualizes the chi-squared values for significant metabolism-cluster associations.
# It uses log transformation for better visualization and includes clustering on both rows and columns.
heatmap_metabolism_clusters <-
  chi2_results_metabolism$chi2_value_matrix[metabolism_clusters_mapping,
                                            names(cluster_mapping_metabolism)]

pheatmap(
  log2(heatmap_metabolism_clusters + 1),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE,
  labels_col = cluster_mapping_metabolism,
  fontsize = 12,             
  fontsize_row = 12,          
  fontsize_col = 12
)

################################################################################
# Visualization of all significant results
# This block is similar to the previous visualization block but is intended for all significant results.
# It prepares the data and generates a heatmap for all significant metabolism-cluster associations.
extract_data(chi2_results_metabolism,  
             tissues_clusters[10:27],
             !rownames(chi2_results_metabolism$p_value_matrix) %in% tissues_clusters) %>%
  rename(cluster = Var1, metabolism = Var2) %>% 
  mutate(fdr = p.adjust(p_value, method = "fdr")) %>% 
  arrange(fdr) %>% 
  filter(fdr < 0.05) %>% 
  filter(number_overlap > 1) %>% 
  .$metabolism %>% 
  as.character() %>% 
  unique() %>% 
  create_mapping_vector() -> metabolism_clusters_mapping

heatmap_metabolism_clusters <-
  chi2_results_metabolism$chi2_value_matrix[metabolism_clusters_mapping,
                                            names(cluster_mapping_metabolism)]

pheatmap(
  log2(heatmap_metabolism_clusters + 1),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE,
  labels_col = cluster_mapping_metabolism
)


################################################################################
##########################################
# 2. tissues vs metabolism visualization #
##########################################
# Extract significant results for tissues vs metabolism
extract_data(
  chi2_results_metabolism,!rownames(chi2_results_metabolism$p_value_matrix) %in% tissues_clusters,
  tissues_clusters[1:9]
) %>%
  rename(cluster = Var1, metabolism = Var2) %>%
  mutate(fdr = p.adjust(p_value, method = "fdr")) %>%
  arrange(fdr) %>%
  filter(fdr < 0.05) %>%
  filter(number_overlap > 1)

# prepare genes overlapping with metabolism
extract_data(
  chi2_results_metabolism,!rownames(chi2_results_metabolism$p_value_matrix) %in% tissues_clusters,
  tissues_clusters[1:9]
) %>%
  rename(cluster = Var1, metabolism = Var2) %>%
  mutate(fdr = p.adjust(p_value, method = "fdr")) %>%
  arrange(fdr) %>%
  filter(fdr < 0.05) %>%
  filter(number_overlap > 1) %>% # chose results where more than one gene was overlapping
  .$overlap_genes %>%
  strsplit(., split = ",") %>%
  unlist %>%
  unique() %>%
  walk(., ~ cat(.x, "\n"))

################################################################################
# prepare row to filter unique gene list
extract_data(chi2_results_metabolism,  
             tissues_clusters[1:9],
             !rownames(chi2_results_metabolism$p_value_matrix) %in% tissues_clusters) %>%
  rename(tissue = Var1, metabolism = Var2) %>% 
  mutate(fdr = p.adjust(p_value, method = "fdr")) %>% 
  arrange(fdr) %>% 
  filter(fdr < 0.05) %>% 
  filter(number_overlap > 1) %>% 
  filter(!(metabolism %in% c("X-17337", "X-11261", "X-11564", "X-21735", "X-11381"))) %>%
  group_by(overlap_genes, tissue) %>% 
  nest() %>% 
  mutate(data = map(data, ~ .x %>% 
                      arrange(fdr) %>% 
                      slice(1))) %>% 
  unnest() %>% 
  ungroup %>% 
  as.data.frame() %>% 
  .$metabolism %>% 
  as.character() %>% 
  unique() %>% 
  create_mapping_vector() -> metabolism_tissue_mapping

heatmap_metabolism_tissue <-
  chi2_results_metabolism$chi2_value_matrix[metabolism_tissue_mapping,
                                            tissues_clusters[1:9]]

pheatmap(
  log2(heatmap_metabolism_tissue + 1),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE,
  fontsize = 12,             
  fontsize_row = 12,          
  fontsize_col = 12
)


# prepare row to filter unique gene list
extract_data(chi2_results_metabolism,  
             tissues_clusters[1:9],
             !rownames(chi2_results_metabolism$p_value_matrix) %in% tissues_clusters) %>%
  rename(cluster = Var1, metabolism = Var2) %>% 
  mutate(fdr = p.adjust(p_value, method = "fdr")) %>% 
  arrange(fdr) %>% 
  filter(fdr < 0.05) %>% 
  filter(number_overlap > 1) %>% 
  .$metabolism %>% 
  as.character() %>% 
  unique() %>% 
  create_mapping_vector() -> metabolism_tissue_mapping

heatmap_metabolism_tissue <-
  chi2_results_metabolism$chi2_value_matrix[metabolism_tissue_mapping,
                                            tissues_clusters[1:9]]

pheatmap(
  log2(heatmap_metabolism_tissue + 1),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE,
)

