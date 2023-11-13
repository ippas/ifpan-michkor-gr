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


##########################################
# 3. visualization for significant results
##########################################
# Extract significant results for clusters vs papers
extract_data(
  chi2_results_phenotypes,
  !rownames(chi2_results_phenotypes$p_value_matrix) %in% tissues_clusters,
  tissues_clusters[10:27]
) %>% 
  dplyr::rename(phenotypes = Var1, cluster = Var2) %>% 
  mutate(fdr = p.adjust(p_value, method = "fdr")) %>%
  arrange(fdr) %>%
  filter(fdr < 0.05) %>% 
  filter(number_overlap > 1) %>% 
  group_by(overlap_genes, cluster) %>% 
  nest() %>% 
  dplyr::mutate(data = map(data, ~ .x %>% 
                      arrange(fdr) %>% 
                      dplyr::slice(1))) %>% 
  unnest() %>% 
  ungroup %>% 
  as.data.frame() %>% 
  .$phenotypes %>% 
  unique %>% 
  as.character() -> significant_phenotypes

# Extract the overlapping genes from the significant results
extract_data(
  chi2_results_phenotypes,
  !rownames(chi2_results_phenotypes$p_value_matrix) %in% tissues_clusters,
  tissues_clusters[10:27]
) %>%
  rename(phenotypes = Var1, cluster = Var2) %>% 
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
  .$overlap_genes %>% 
  strsplit(., split = ",") %>% 
  unlist %>% 
  unique() %>% 
  walk(., ~cat(.x, "\n"))

# prepare data to heatmap
heatmap_phenotypes_clusters <-
  chi2_results_phenotypes$chi2_value_matrix[significant_phenotypes,
                                            tissues_clusters[10:27]]

# create heatmap
pheatmap(
  log2(heatmap_phenotypes_clusters + 1),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE,
  fontsize = 12,             
  fontsize_row = 15,          
  fontsize_col = 15
)

###############################################################################
# tissue
extract_data(
  chi2_results_phenotypes,
  !rownames(chi2_results_phenotypes$p_value_matrix) %in% tissues_clusters,
  tissues_clusters[1:9]
) %>%
  rename(phenotypes = Var1, cluster = Var2) %>% 
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
  .$phenotypes %>% 
  unique %>% 
  as.character() -> significant_phenotypes

# prepare data to heatmap
heatmap_phenotypes_tissues <-
  chi2_results_phenotypes$chi2_value_matrix[significant_phenotypes,
                                            tissues_clusters[1:9]]

# create heatmap
pheatmap(
  log2(heatmap_phenotypes_tissues + 1),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE,
  fontsize = 12,             
  fontsize_row = 15,          
  fontsize_col = 15
)

#################################################################
# 4. visualization for significant results with additional filter 
#################################################################
phenotypes_biobank_genes_list %>%
  # Keep only those vectors that have an intersection of at least two genes with master_genes
  keep(~ length(intersect(.x, {genes_list %>% unlist() %>% unique})) > 1) -> phenotypes_filt_genes_list



chi2_results_phenotypes_filt <-
  perform_chi2_tests(c(papers_gene_list[tissues_clusters], phenotypes_filt_genes_list),
                     hgnc_symbols_vector_v110)

extract_data(
  chi2_results_phenotypes_filt,
  !rownames(chi2_results_phenotypes_filt$p_value_matrix) %in% tissues_clusters,
  tissues_clusters[10:27]
) %>%
  rename(phenotypes = Var1, cluster = Var2) %>% 
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
  .$phenotypes %>% 
  as.character %>%
  unique -> significant_phenotypes

# Extract the overlapping genes from the significant results
extract_data(
  chi2_results_phenotypes_filt,
  !rownames(chi2_results_phenotypes_filt$p_value_matrix) %in% tissues_clusters,
  tissues_clusters[10:27]
) %>%
  rename(phenotypes = Var1, cluster = Var2) %>% 
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
  .$overlap_genes %>% 
  strsplit(., split = ",") %>% 
  unlist %>% 
  unique() %>% 
  walk(., ~cat(.x, "\n"))

# prepare data to heatmap
heatmap_phenotypes_clusters <-
  chi2_results_phenotypes_filt$chi2_value_matrix[significant_phenotypes,
                                                 tissues_clusters[10:27]]

# create heatmap
pheatmap(
  log2(heatmap_phenotypes_clusters + 1),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE,
  fontsize = 12,             
  fontsize_row = 15,          
  fontsize_col = 15
)

################################################################################
# tissue
extract_data(
  chi2_results_phenotypes_filt,
  !rownames(chi2_results_phenotypes_filt$p_value_matrix) %in% tissues_clusters,
  tissues_clusters[1:9]
) %>%
  rename(phenotypes = Var1, cluster = Var2) %>% 
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
  .$phenotypes %>% 
  as.character %>%
  unique -> significant_phenotypes


# prepare data to heatmap
heatmap_phenotypes_tissues <-
  chi2_results_phenotypes_filt$chi2_value_matrix[significant_phenotypes,
                                                 tissues_clusters[1:9]]

# create heatmap
pheatmap(
  log2(heatmap_phenotypes_tissues + 1),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE,
  fontsize = 12,             
  fontsize_row = 15,          
  fontsize_col = 15
)
