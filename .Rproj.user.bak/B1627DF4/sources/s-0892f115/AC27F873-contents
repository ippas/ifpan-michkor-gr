# function
create_mapping_vector <- function(input_vector) {
  # Extract the last element after the underscore
  extracted_elements <- sapply(strsplit(input_vector, "_"), tail, 1)
  
  # Format the extracted elements
  # formatted_elements <- paste("cluster", extracted_elements)
  
  # Create the named vector
  mapping_vector <- setNames(extracted_elements, input_vector)
  
  return(mapping_vector)
}

#####################################
# 1. Cluster vs metabolism analysis #
#####################################
# prepare data to chi2
metabolism_preprocesing_data <- gr_gene_database %>%
  filter(
    grepl("omicspred_metabolon", source) |
      grepl("omicspred_nithingale_", source) |
      grepl("marpiech", gene_list_index)
  ) %>%
  filter(gene_list_index != "all_significant_genes_marpiech") %>%
  mutate(source = ifelse(source == "pmid:NA", "marpiech_tissues_dex", source)) %>%
  extract_keys_values("info", c("BiomarkerName", "BiochemicalName")) %>%
  mutate(label = ifelse(is.na(BiochemicalName), ifelse(is.na(BiomarkerName), gene_list_index, BiomarkerName), BiochemicalName)) %>%
  mutate(label = ifelse(source == "marpiech_tissues_dex", label, paste0(gene_list_index, "_", label))) %>%
  select(c(source, gene_list_index, gene_list_number, hgnc_symbol, label)) %>%
  filter(!is.na(hgnc_symbol)) %>% 
  distinct() %>%
  filter(!(hgnc_symbol %in% hgnc_to_remove)) %>% 
  group_by(gene_list_number) %>%
  nest() %>%
  mutate(n_genes = map_int(data, ~ nrow(.))) %>% 
  unnest(cols = c(data)) %>% 
  ungroup() %>% 
  filter(n_genes > 2)

# check data
metabolism_preprocesing_data %>% tail

# prepare list of genes for metabolism genes
split(metabolism_preprocesing_data$hgnc_symbol, metabolism_preprocesing_data$label) -> metabolism_gene_list

# calculate chi2
results_metabolism <- perform_chi2_tests(metabolism_gene_list, hgnc_symbols_vector)

# filter significant results
# fdr threshold = 0.05


# prepare genes overlapping with metabolomics
extract_data(results_metabolism,  1:18, 19:ncol(results_metabolism$p_value_matrix)) %>%
  rename(cluster = Var1, metabolism = Var2) %>% 
  mutate(fdr = p.adjust(p_value, method = "fdr")) %>% 
  arrange(fdr) %>% 
  filter(fdr < 0.05) %>% 
  filter(number_overlap > 1) %>% # chose results where more than one gene was overlapping
  .$overlap_genes %>% 
  strsplit(., split = ",") %>% 
  unlist %>% 
  unique() %>% 
  walk(., ~cat(.x, "\n"))

# filter significant results
extract_data(results_metabolism,  1:18, 19:ncol(results_metabolism$p_value_matrix)) %>%
  rename(cluster = Var1, metabolism = Var2) %>% 
  mutate(fdr = p.adjust(p_value, method = "fdr")) %>% 
  arrange(fdr) %>% 
  filter(fdr < 0.05) %>% 
  filter(number_overlap > 1)

# prepare data to heatmap visualisation 
# prepare column to filter 
cluster_mapping_metabolism <- c("marpiech_9" = "cluster 9", 
                                 "marpiech_12" = "cluster 12", 
                                 "marpiech_tissues_dex_up" = "upregulated genes")

# prepare row to filter
extract_data(results_metabolism,  1:18, 19:ncol(results_metabolism$p_value_matrix)) %>%
  rename(cluster = Var1, metabolism = Var2) %>% 
  mutate(fdr = p.adjust(p_value, method = "fdr")) %>% 
  arrange(fdr) %>% 
  filter(fdr < 0.05) %>% 
  filter(number_overlap > 1) %>% 
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
  create_mapping_vector() -> metabolism_to_heatmap

results_metabolism$chi2_value_matrix[names(metabolism_to_heatmap), names(cluster_mapping_metabolism)] -> metabolism_data_heatmap

pheatmap(
  log2(metabolism_data_heatmap + 1),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE,
  labels_col = cluster_mapping_metabolism[colnames(metabolism_data_heatmap)],
  labels_row = metabolism_to_heatmap,
  fontsize = 12,             
  fontsize_row = 12,          
  fontsize_col = 12
)


#####################################
# 2. Tissue vs metabolism analysis #
#####################################
# grepl("omicspred_metabolon", source) |
#   grepl("omicspred_nithingale_", source) |
metabolism_tissue_preprocesing_data <- gr_gene_database %>%
  filter(
    grepl("omicspred_nithingale_", source) |
      grepl("omicspred_metabolon", source)  |
      gene_list_index == "all_significant_genes_marpiech"
  ) %>%
  mutate(source = ifelse(source == "pmid:NA", "marpiech_tissues_dex", source)) %>%
  extract_keys_values("info", c("BiomarkerName", "BiochemicalName")) %>% 
  mutate(label = ifelse(!is.na(BiochemicalName), BiochemicalName, ifelse(!is.na(BiomarkerName), BiomarkerName, ifelse(!is.na(tissue), tissue, NA)))) %>% 
  mutate(label = ifelse(source == "marpiech_tissues_dex", label, paste0(gene_list_index, "_", label))) %>%
  select(c(source, gene_list_index, gene_list_number, hgnc_symbol, label)) %>%
  filter(!is.na(hgnc_symbol)) %>% 
  distinct() %>%
  filter(!(hgnc_symbol %in% hgnc_to_remove)) %>% 
  group_by(gene_list_number) %>%
  nest() %>%
  mutate(n_genes = map_int(data, ~ nrow(.))) %>% 
  unnest(cols = c(data)) %>% 
  ungroup() %>% 
  filter(n_genes > 2)

# prepare list of genes for metabolism genes
split(metabolism_tissue_preprocesing_data$hgnc_symbol, metabolism_tissue_preprocesing_data$label) -> metabolism_tissue_gene_list

# calculate chi2
results_tissue_metabolism <- perform_chi2_tests(metabolism_tissue_gene_list, hgnc_symbols_vector)

extract_data(results_tissue_metabolism,  1:18, 19:ncol(results_tissue_metabolism$p_value_matrix)) %>%
  rename(cluster = Var1, metabolism = Var2) %>% 
  mutate(fdr = p.adjust(p_value, method = "fdr")) %>% 
  arrange(fdr) %>% 
  filter(fdr < 0.05) %>% 
  filter(number_overlap > 1)


# testing code
tissues <-c(
  "adrenal-cortex",
  "perigonadal-adipose-tissue",
  "hypothalamus",
  "kidneys",
  "liver",
  "anterior-thigh",
  "pituitary-gland",
  "spleen",
  "lung"
)

# filter significant results
extract_data(results_tissue_metabolism, tissues, setdiff(colnames(results_tissue_metabolism$chi2_value_matrix), tissues)) %>% 
  rename(tissue = Var1, metabolism = Var2) %>% 
  mutate(fdr = p.adjust(p_value, method = "fdr")) %>% 
  arrange(fdr) %>% 
  filter(fdr < 0.05) %>% 
  filter(number_overlap > 1)


# prepare genes overlapping with metabolomics
extract_data(results_tissue_metabolism, tissues, setdiff(colnames(results_tissue_metabolism$chi2_value_matrix), tissues)) %>% 
  rename(tissue = Var1, metabolism = Var2) %>% 
  mutate(fdr = p.adjust(p_value, method = "fdr")) %>% 
  arrange(fdr) %>% 
  filter(fdr < 0.05) %>% 
  filter(number_overlap > 1) %>%   
  .$overlap_genes %>% 
  strsplit(., split = ",") %>% 
  unlist %>% 
  unique() %>% 
  walk(., ~cat(.x, "\n"))

# prepare row to filter
extract_data(results_tissue_metabolism, tissues, setdiff(colnames(results_tissue_metabolism$chi2_value_matrix), tissues)) %>% 
  rename(tissue = Var1, metabolism = Var2) %>% 
  mutate(fdr = p.adjust(p_value, method = "fdr")) %>% 
  arrange(fdr) %>% 
  filter(fdr < 0.05) %>% 
  filter(number_overlap > 1) %>% 
  group_by(overlap_genes) %>% 
  nest() %>% 
  mutate(data = map(data, ~ .x %>% 
                      arrange(fdr) %>% 
                      slice(1))) %>% 
  unnest() %>% 
  ungroup %>% 
  as.data.frame() %>% 
  .$metabolism %>% 
  as.character() %>% unique() %>% 
  create_mapping_vector() -> metabolism_tissue_to_heatmap



results_tissue_metabolism$chi2_value_matrix[names(metabolism_tissue_to_heatmap), tissues]  -> metabolism_tissue_data_heatmap

pheatmap(
  log2(metabolism_tissue_data_heatmap + 1),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE,
  labels_row = metabolism_tissue_to_heatmap,
  fontsize = 12,             
  fontsize_row = 12,          
  fontsize_col = 12
)

## testing
metabolism_tissue_gene_list$Metabolon_OPGS003129_model_X-11564

gr_gene_database %>%
  filter(gene_list_index == "Metabolon_OPGS003129_model") %>% .$gene_name %>% print(quote = F) 
  # tail(4000) %>%  head

# problem z listatmi, bo z jednej strony te same geny, ale później pojawia się
# problem bo mogą wylecieć jakieś inne nazwy procesów, 
# oraz czasem są dziwne nazwy (numeryczne), przez co nie wiadomo co to jest
# wyświetlanie wszystkich też nie do końca dobre bo będzie jeszcze mniej czytelnie na heatmapie
extract_data(results_tissue_metabolism, tissues, setdiff(colnames(results_tissue_metabolism$chi2_value_matrix), tissues)) %>% 
  rename(tissue = Var1, metabolism = Var2) %>% 
  mutate(fdr = p.adjust(p_value, method = "fdr")) %>% 
  arrange(fdr) %>% 
  filter(fdr < 0.05) %>% 
  filter(number_overlap > 1) %>% 
  filter(overlap_genes == "SLC16A9,IRF1") %>% 
  select(c(metabolism, overlap_genes)) %>% unique()
