# refine_and_label_gene_lists(papers_data_preprocessing, columns = c("source", "simple_tissue", "cell", "environment", "treatment", "dose", "time")) -> tmp 
secondary_gene_lists_all <- list()

refine_and_label_gene_lists(papers_data_preprocessing, columns = c("source")) -> secondary_gene_lists_all$refine_gene_lists

summarize_genes_distribution(secondary_gene_lists_all$refine_gene_lists) -> secondary_gene_lists_all$genes_distribution

# Example of cumulative sum threshold values
cumsum_thresholds <- c(0.99, 0.995, 0.996, 0.997, 0.998, 0.999)

# Calling the function
secondary_gene_lists_all$master_gene_lists <- create_master_gene_lists(cumsum_thresholds, secondary_gene_lists_all$genes_distribution$genes_frequency)

# Example of cumulative sum threshold values
freq_thresholds <- list(1, 2, 3, c(1, 2), c(1, 2, 3))

secondary_gene_lists_all$rare_gene_lists <- create_rare_gene_lists(freq_thresholds = freq_thresholds, 
                                          gene_distribution = secondary_gene_lists_all$genes_distribution$genes_frequency)



papers_data_preprocessing %>% 
  filter(!(treatment %in% c("TNF", "LPS", "vehicle-ethanol", "TNFalpha"))) %>% 
  filter(!(source %in% c("marpiech_tissues_dex",  "pmid:37217509"))) %>% 
  filter(source != "marpiech_tissues") %>% 
  mutate(label = paste(source, tissue, cell, treatment, treatment_type, regulation, sep = "_")) %>% 
  select(c(hgnc_symbol, label)) -> tmp_df

tmp_df %>% 
  # filter((hgnc_symbol %in% secondary_gene_lists_all$rare_gene_lists$freq_123)) %>%
  split(.$label, .$hgnc_symbol) %>% 
  lapply(., unique) %>% 
  lapply(., function(x){x$hgnc_symbol}) %>% 
  lapply(., length) %>% as.data.frame() %>% t %>% set_colnames("n_origin") %>% as.data.frame() %>% 
  rownames_to_column(var = "list_name") -> origin_size


tmp_df %>% 
  filter(!(hgnc_symbol %in% secondary_gene_lists_all$rare_gene_lists$freq_1)) %>%
  split(.$label, .$hgnc_symbol) %>% 
  lapply(., unique) %>% 
  lapply(., function(x){x$hgnc_symbol}) %>% 
  lapply(., length) %>% as.data.frame() %>% t %>% set_colnames("n_minus1") %>% as.data.frame() %>% 
  rownames_to_column(var = "list_name") -> minus1_size

tmp_df %>% 
  filter(!(hgnc_symbol %in% secondary_gene_lists_all$rare_gene_lists$freq_12)) %>%
  split(.$label, .$hgnc_symbol) %>% 
  lapply(., unique) %>% 
  lapply(., function(x){x$hgnc_symbol}) %>% 
  lapply(., length) %>% as.data.frame() %>% t %>% set_colnames("n_minus12") %>% as.data.frame() %>% 
  rownames_to_column(var = "list_name") -> minus12_size

tmp_df %>% 
  filter(!(hgnc_symbol %in% secondary_gene_lists_all$rare_gene_lists$freq_123)) %>%
  split(.$label, .$hgnc_symbol) %>% 
  lapply(., unique) %>% 
  lapply(., function(x){x$hgnc_symbol}) %>% 
  lapply(., length) %>% as.data.frame() %>% t %>% set_colnames("n_minus123") %>% as.data.frame() %>% 
  rownames_to_column(var = "list_name") -> minus123_size

left_join(origin_size, minus1_size, by = "list_name") %>% 
  left_join(., minus12_size, by = "list_name") %>% 
  left_join(., minus123_size, by = "list_name") %>% 
  gather(key = "group", value = "n_genes", -list_name) %>% 
  mutate(group = factor(group, levels = c("n_origin", "n_minus1", "n_minus12", "n_minus123"))) %>% 
  ggplot(aes(x = n_genes)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black") +  # adjust the number of bins as needed
  facet_grid(rows = vars(group)) +
  theme_minimal() +
  labs(x = "Number of Genes", y = "Frequency", title = "Histograms of Gene Counts by Group")

# Merge the data frames using left joins and gather to a long format
origin_size %>%
  left_join(minus1_size, by = "list_name") %>%
  left_join(minus12_size, by = "list_name") %>%
  left_join(minus123_size, by = "list_name") %>%
  gather(key = "group", value = "n_genes", -list_name) %>% 
  mutate(group = factor(group, levels = c("n_origin", "n_minus1", "n_minus12", "n_minus123"))) %>% 
  ggplot(aes(x = group, y = n_genes, fill = group)) +
  geom_boxplot() +
  geom_jitter() +
  theme_minimal() +
  labs(x = "Group", y = "Number of Genes", title = "Boxplot of Gene Counts by Group") +
  scale_fill_brewer(palette = "Pastel1")

################################################################################
gr_database_blocked_gene_lists$GR_dependent_gene_lists_pmid_tissue_cell_treatment_type_regulation_minus1$metadata <- 
  papers_data_preprocessing %>% 
  filter(!(treatment %in% c("TNF", "LPS", "vehicle-ethanol", "TNFalpha"))) %>% 
  filter(!(source %in% c("marpiech_tissues_dex",  "pmid:37217509"))) %>% 
  filter(source != "marpiech_tissues") %>% 
  mutate(label = paste(source, tissue, cell, treatment, treatment_type, regulation, sep = "_")) %>%
  filter(!(hgnc_symbol %in% secondary_gene_lists_all$rare_gene_lists$freq_1)) %>% 
  select(-c(index, gene_name, ensembl_gene_id, ensembl_transcript_id, refseq_mrna_id, hgnc_symbol, alias, info, fdr, log2ratio))

gr_database_blocked_gene_lists$GR_dependent_gene_lists_pmid_tissue_cell_treatment_type_regulation_minus1$gene_lists <- papers_data_preprocessing %>% 
  filter(!(treatment %in% c("TNF", "LPS", "vehicle-ethanol", "TNFalpha"))) %>% 
  filter(!(source %in% c("marpiech_tissues_dex",  "pmid:37217509"))) %>% 
  filter(source != "marpiech_tissues") %>% 
  mutate(label = paste(source, tissue, cell, treatment, treatment_type, regulation, sep = "_")) %>% 
  filter(!(hgnc_symbol %in% secondary_gene_lists_all$rare_gene_lists$freq_1)) %>% 
  select(c(hgnc_symbol, label)) %>% 
  split(.$label, .$hgnc_symbol) %>% 
  lapply(., unique) %>% 
  lapply(., function(x){x$hgnc_symbol})


################################################################################
gr_database_blocked_gene_lists$GR_dependent_gene_lists_pmid_tissue_cell_treatment_type_regulation_minus123$metadata <- 
  papers_data_preprocessing %>% 
  filter(!(treatment %in% c("TNF", "LPS", "vehicle-ethanol", "TNFalpha"))) %>% 
  filter(!(source %in% c("marpiech_tissues_dex",  "pmid:37217509"))) %>% 
  filter(source != "marpiech_tissues") %>% 
  mutate(label = paste(source, tissue, cell, treatment, treatment_type, regulation, sep = "_")) %>%
  filter(!(hgnc_symbol %in% secondary_gene_lists_all$rare_gene_lists$freq_123)) %>% 
  select(-c(index, gene_name, ensembl_gene_id, ensembl_transcript_id, refseq_mrna_id, hgnc_symbol, alias, info, fdr, log2ratio))

gr_database_blocked_gene_lists$GR_dependent_gene_lists_pmid_tissue_cell_treatment_type_regulation_minus123$gene_lists <- papers_data_preprocessing %>% 
  filter(!(treatment %in% c("TNF", "LPS", "vehicle-ethanol", "TNFalpha"))) %>% 
  filter(!(source %in% c("marpiech_tissues_dex",  "pmid:37217509"))) %>% 
  filter(source != "marpiech_tissues") %>% 
  mutate(label = paste(source, tissue, cell, treatment, treatment_type, regulation, sep = "_")) %>% 
  filter(!(hgnc_symbol %in% secondary_gene_lists_all$rare_gene_lists$freq_123)) %>% 
  select(c(hgnc_symbol, label)) %>% 
  split(.$label, .$hgnc_symbol) %>% 
  lapply(., unique) %>% 
  lapply(., function(x){x$hgnc_symbol})


################################################################################
# mclusters vs papers
analyze_gene_list_overlap(
  row_lists = gr_database_blocked_gene_lists$marpiech_cluster_dex,
  col_lists = gr_database_blocked_gene_lists$GR_dependent_gene_lists_pmid_tissue_cell_treatment_type_regulation_minus1$gene_lists,
  reference_hgnc_vector = hgnc_symbols_vector_v110, 
  keep_original_data = TRUE, 
  fdr_threshold = 0.01, 
  overlap_threshold = 3
) -> restore_overlap_analysis$mcluster_vs_papers_minus1

# gr_database_blocked_gene_lists$GR_dependent_gene_lists_pmid_tissue_cell_treatment_type_regulation$gene_lists %>%
#   lapply(., length) %>%
#   unname() %>%
#   unlist() %>%
#   hist(xlab = "Length of Gene Lists", ylab = "Frequency", main = "Histogram of Gene List Lengths")

restore_overlap_analysis$mcluster_vs_papers_minus1$significant_uniq_data$overlap_genes

fast_heatmap(data = restore_overlap_analysis$mcluster_vs_papers_minus1,
             data_type = "significant_uniq_data")

################################################################################
# mclusters vs papers
analyze_gene_list_overlap(
  row_lists = gr_database_blocked_gene_lists$marpiech_cluster_dex,
  col_lists = gr_database_blocked_gene_lists$GR_dependent_gene_lists_pmid_tissue_cell_treatment_type_regulation_minus123$gene_lists,
  reference_hgnc_vector = hgnc_symbols_vector_v110, 
  keep_original_data = TRUE, 
  fdr_threshold = 0.01, 
  overlap_threshold = 3
) -> restore_overlap_analysis$mcluster_vs_papers_minus123



fast_heatmap(data = restore_overlap_analysis$mcluster_vs_papers_minus123,
             data_type = "significant_uniq_data")

################################################################################
gr_database_blocked_gene_lists$GR_dependent_gene_lists_pmid_tissue_cell_treatment_type_regulation_only1$metadata <- 
  papers_data_preprocessing %>% 
  filter(!(treatment %in% c("TNF", "LPS", "vehicle-ethanol", "TNFalpha"))) %>% 
  filter(!(source %in% c("marpiech_tissues_dex",  "pmid:37217509"))) %>% 
  filter(source != "marpiech_tissues") %>% 
  mutate(label = paste(source, tissue, cell, treatment, treatment_type, regulation, sep = "_")) %>%
  filter((hgnc_symbol %in% secondary_gene_lists_all$rare_gene_lists$freq_1)) %>% 
  select(-c(index, gene_name, ensembl_gene_id, ensembl_transcript_id, refseq_mrna_id, hgnc_symbol, alias, info, fdr, log2ratio))

gr_database_blocked_gene_lists$GR_dependent_gene_lists_pmid_tissue_cell_treatment_type_regulation_only1$gene_lists <- papers_data_preprocessing %>% 
  filter(!(treatment %in% c("TNF", "LPS", "vehicle-ethanol", "TNFalpha"))) %>% 
  filter(!(source %in% c("marpiech_tissues_dex",  "pmid:37217509"))) %>% 
  filter(source != "marpiech_tissues") %>% 
  mutate(label = paste(source, tissue, cell, treatment, treatment_type, regulation, sep = "_")) %>% 
  filter((hgnc_symbol %in% secondary_gene_lists_all$rare_gene_lists$freq_1)) %>% 
  select(c(hgnc_symbol, label)) %>% 
  split(.$label, .$hgnc_symbol) %>% 
  lapply(., unique) %>% 
  lapply(., function(x){x$hgnc_symbol})


################################################################################
# mclusters vs papers
analyze_gene_list_overlap(
  row_lists = gr_database_blocked_gene_lists$marpiech_cluster_dex,
  col_lists = gr_database_blocked_gene_lists$GR_dependent_gene_lists_pmid_tissue_cell_treatment_type_regulation_only1$gene_lists,
  reference_hgnc_vector = hgnc_symbols_vector_v110, 
  keep_original_data = TRUE, 
  fdr_threshold = 0.01, 
  overlap_threshold = 0
) -> restore_overlap_analysis$mcluster_vs_papers_only1

fast_heatmap(data = restore_overlap_analysis$mcluster_vs_papers_minus1,
             data_type = "significant_uniq_data")


################################################################################
gr_database_blocked_gene_lists$GR_dependent_gene_lists_pmid_tissue_cell_treatment_type_regulation_only123$metadata <- 
  papers_data_preprocessing %>% 
  filter(!(treatment %in% c("TNF", "LPS", "vehicle-ethanol", "TNFalpha"))) %>% 
  filter(!(source %in% c("marpiech_tissues_dex",  "pmid:37217509"))) %>% 
  filter(source != "marpiech_tissues") %>% 
  mutate(label = paste(source, tissue, cell, treatment, treatment_type, regulation, sep = "_")) %>%
  filter((hgnc_symbol %in% secondary_gene_lists_all$rare_gene_lists$freq_123)) %>% 
  select(-c(index, gene_name, ensembl_gene_id, ensembl_transcript_id, refseq_mrna_id, hgnc_symbol, alias, info, fdr, log2ratio))

gr_database_blocked_gene_lists$GR_dependent_gene_lists_pmid_tissue_cell_treatment_type_regulation_only123$gene_lists <- papers_data_preprocessing %>% 
  filter(!(treatment %in% c("TNF", "LPS", "vehicle-ethanol", "TNFalpha"))) %>% 
  filter(!(source %in% c("marpiech_tissues_dex",  "pmid:37217509"))) %>% 
  filter(source != "marpiech_tissues") %>% 
  mutate(label = paste(source, tissue, cell, treatment, treatment_type, regulation, sep = "_")) %>% 
  filter((hgnc_symbol %in% secondary_gene_lists_all$rare_gene_lists$freq_123)) %>% 
  select(c(hgnc_symbol, label)) %>% 
  split(.$label, .$hgnc_symbol) %>% 
  lapply(., unique) %>% 
  lapply(., function(x){x$hgnc_symbol})


# mclusters vs papers
analyze_gene_list_overlap(
  row_lists = gr_database_blocked_gene_lists$marpiech_cluster_dex,
  col_lists = gr_database_blocked_gene_lists$GR_dependent_gene_lists_pmid_tissue_cell_treatment_type_regulation_only123$gene_lists,
  reference_hgnc_vector = hgnc_symbols_vector_v110, 
  keep_original_data = TRUE, 
  fdr_threshold = 0.1, 
  overlap_threshold = 2
) -> restore_overlap_analysis$mcluster_vs_papers_only123

fast_heatmap(data = restore_overlap_analysis$mcluster_vs_papers_only123,
             data_type = "significant_uniq_data")
