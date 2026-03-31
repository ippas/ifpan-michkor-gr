

set_data$gr_source$rare_gene_lists$freq_123 %>% length()
lapply(set_data$gr_source$rare_gene_lists, length)


gr_database_blocked_gene_lists$marpiech_cluster_dex %>% 
  unname() %>% unlist %>% unique %>% 
  intersect(., set_data$gr_source$rare_gene_lists$freq_123)


papers_data_preprocessing %>% 
  filter(!(treatment %in% treatment_to_remove)) %>% 
  filter(!(hgnc_symbol %in% set_data$gr_source$rare_gene_lists$freq_123)) %>% 
  refine_gene_lists(data = .,
                  columns = c("source", "tissue", "cell", "treatment", "treatment_type", "regulation", "comparison", "environment"),
                  cumsum_thresholds = cumsum_thresholds,
                  freq_thresholds = freq_thresholds,
                  keep_column = c("system", "simple_tissue", "method")) -> set_data$gr_minus123

papers_data_preprocessing %>% 
  filter(!(treatment %in% treatment_to_remove)) %>% 
  filter((hgnc_symbol %in% set_data$gr_source$rare_gene_lists$freq_123)) %>% 
  refine_gene_lists(data = .,
                    columns = c("source", "tissue", "cell", "treatment", "treatment_type", "regulation", "comparison", "environment"),
                    cumsum_thresholds = cumsum_thresholds,
                    freq_thresholds = freq_thresholds,
                    keep_column = c("system", "simple_tissue", "method")) -> set_data$gr_only123


papers_data_preprocessing %>% 
  filter(!(treatment %in% treatment_to_remove)) %>% 
  filter((hgnc_symbol %in% set_data$gr_source$rare_gene_lists$freq_123)) %>% 
  refine_gene_lists(data = .,
                    columns = c("source"),
                    cumsum_thresholds = cumsum_thresholds,
                    freq_thresholds = freq_thresholds) -> set_data$gr_source_only123 



################################################################################
# investigation of size gene lists
set_data$gr_minus123$refine_gene_lists$metadata %>% dim
  
set_data$gr_minus123$refine_gene_lists$metadata %>% 
  ggplot(aes(x = gene_count)) + 
  geom_histogram() +
  theme_minimal()


set_data$gr_minus123$refine_gene_lists$metadata %>% 
  ggplot(aes(x = gene_count)) + 
  geom_histogram() +
  facet_wrap(~ method) +
  theme_minimal()

set_data$gr_minus123$refine_gene_lists$metadata %>% 
  filter(regulation %in% c("up", "down")) %>% 
  ggplot(aes(x = gene_count)) + 
  geom_histogram() +
  facet_wrap(~ regulation) +
  theme_minimal()

set_data$gr_minus123$refine_gene_lists$metadata %>% 
  mutate(decile = ntile(gene_count, 10)) %>%
  group_by(decile) %>%
  summarize(
    mean_gene_count = mean(gene_count),
    median_gene_count = median(gene_count),
    min_gene_count = min(gene_count),
    max_gene_count = max(gene_count),
    sd_gene_count = sd(gene_count)
  )  

set_data$gr_minus123$refine_gene_lists$metadata$gene_count %>% summary


################################################################################
# overlap analysis: clusters vs papers
analyze_gene_list_overlap(
  row_lists = gr_database_blocked_gene_lists$marpiech_cluster_dex,
  col_lists = set_data$gr_minus123$refine_gene_lists$gene_lists,
  reference_hgnc_vector = hgnc_symbols_vector_v110, 
  keep_original_data = TRUE, 
  fdr_threshold = 0.01, 
  overlap_threshold = 3
) -> results_overlap$gr_minus123$mclusters

fast_heatmap(data =   results_overlap$gr_minus123$mclusters,
             data_type = "significant_uniq_data")

# overlap analysis: clusters vs papers
analyze_gene_list_overlap(
  row_lists = filter_list_by_names( gr_database_blocked_gene_lists$marpiech_cluster_dex, c("cluster_18", "cluster_17"), operation = "!%in%"),
  col_lists = set_data$gr_minus123$refine_gene_lists$gene_lists,
  reference_hgnc_vector = hgnc_symbols_vector_v110, 
  keep_original_data = TRUE, 
  fdr_threshold = 0.01, 
  overlap_threshold = 3
) -> results_overlap$gr_minus123$mclusters_minusUpDown

fast_heatmap(data = results_overlap$gr_minus123$mclusters_minusUpDown,
             data_type = "significant_uniq_data")


# overlap analysis: clusters vs papers
analyze_gene_list_overlap(
  row_lists = filter_list_by_names( gr_database_blocked_gene_lists$marpiech_cluster_dex, c("cluster_18", "cluster_17"), operation = "!%in%"),
  col_lists = set_data$gr_minus123$refine_gene_lists$gene_lists,
  reference_hgnc_vector = hgnc_symbols_vector_v110, 
  keep_original_data = TRUE, 
  fdr_threshold = 0.01, 
  overlap_threshold = 10
) -> results_overlap$gr_minus123$mclusters_minusUpDown_overlap10

fast_heatmap(data = results_overlap$gr_minus123$mclusters_minusUpDown_overlap10,
             data_type = "significant_uniq_data", overlap_threshold = 10)

################################################################################


set_data$gr_source$refine_gene_lists$metadata





################################################################################
# overlap analysis: clusters vs papers
analyze_gene_list_overlap(
  row_lists = gr_database_blocked_gene_lists$marpiech_cluster_dex,
  col_lists = set_data$gr_minus123$refine_gene_lists$gene_lists,
  reference_hgnc_vector = hgnc_symbols_vector_v110, 
  keep_original_data = TRUE, 
  fdr_threshold = 0.01, 
  overlap_threshold = 3
) -> results_overlap$gr_minus123$mclusters

fast_heatmap(data = results_overlap$gr_minus123$mclusters,
             data_type = "significant_uniq_data")

summary_overlap_results(results_overlap$gr_minus123$mclusters$significant_uniq_data)

# overlap analysis: clusters (without UP and DOWN) vs papers
analyze_gene_list_overlap(
  row_lists = filter_list_by_names( gr_database_blocked_gene_lists$marpiech_cluster_dex, c("cluster_18", "cluster_17"), operation = "!%in%"),
  col_lists = set_data$gr_minus123$refine_gene_lists$gene_lists,
  reference_hgnc_vector = hgnc_symbols_vector_v110, 
  keep_original_data = TRUE, 
  fdr_threshold = 0.01, 
  overlap_threshold = 3
) -> results_overlap$gr_minus123$mclusters_minusUpDown

fast_heatmap(data =   results_overlap$gr_minus123$mclusters_minusUpDown,
             data_type = "significant_uniq_data", overlap_threshold = 3)

summary_overlap_results(results_overlap$gr_minus123$mclusters_minusUpDown$significant_uniq_data)

# overlap analysis: clusters (without UP and DOWN) vs papers. (n overlap threshold = 10 genes)
analyze_gene_list_overlap(
  row_lists = filter_list_by_names( gr_database_blocked_gene_lists$marpiech_cluster_dex, c("cluster_18", "cluster_17"), operation = "!%in%"),
  col_lists = set_data$gr_minus123$refine_gene_lists$gene_lists,
  reference_hgnc_vector = hgnc_symbols_vector_v110, 
  keep_original_data = TRUE, 
  fdr_threshold = 0.01, 
  overlap_threshold = 10
) -> results_overlap$gr_minus123$mclusters_minusUpDown_overlap10

fast_heatmap(data =   results_overlap$gr_minus123$mclusters_minusUpDown_overlap10,
             data_type = "significant_uniq_data", overlap_threshold = 10)

summary_overlap_results(results_overlap$gr_minus123$mclusters_minusUpDown_overlap10$significant_uniq_data)



################################################################################
# analysis rare genes
set_data$gr_source_only123$refine_gene_lists$metadata %>% 
  ggplot(aes(x = gene_count)) + 
  geom_histogram() +
  theme_minimal()


# overlap analysis: clusters vs papers
analyze_gene_list_overlap(
  row_lists = gr_database_blocked_gene_lists$marpiech_cluster_dex,
  col_lists =  set_data$gr_only123$refine_gene_lists$gene_lists,
  reference_hgnc_vector = hgnc_symbols_vector_v110, 
  keep_original_data = TRUE, 
  fdr_threshold = 0.01, 
  overlap_threshold = 3
) -> results_overlap$gr_only123$mclusters

fast_heatmap(data =   results_overlap$gr_only123$mclusters,
             data_type = "significant_uniq_data")

results_overlap$gr_only123$mclusters$significant_uniq_data$overlap_genes

results_overlap$gr_only123$mclusters$significant_uniq_data$df


################################################################################
# minus123 and minusCumsum0.995

papers_data_preprocessing %>% 
  filter(!(treatment %in% treatment_to_remove)) %>% 
  filter(!(hgnc_symbol %in% set_data$gr_source$rare_gene_lists$freq_123)) %>% 
  filter(!(hgnc_symbol %in% set_data$gr_source$master_gene_lists$master_cumsum_0.995)) %>% 
  refine_gene_lists(data = .,
                    columns = c("source", "tissue", "cell", "treatment", "treatment_type", "regulation", "comparison", "environment"),
                    cumsum_thresholds = cumsum_thresholds,
                    freq_thresholds = freq_thresholds,
                    keep_column = c("system", "simple_tissue", "method")) -> set_data$gr_minus123_minusCumsum0.995


# overlap analysis: clusters vs papers
analyze_gene_list_overlap(
  row_lists = gr_database_blocked_gene_lists$marpiech_cluster_dex,
  col_lists =  set_data$gr_minus123_minusCumsum0.995$refine_gene_lists$gene_lists,
  reference_hgnc_vector = hgnc_symbols_vector_v110, 
  keep_original_data = TRUE, 
  fdr_threshold = 0.01, 
  overlap_threshold = 3
) -> results_overlap$gr_minus123_minusCumsum0.995$mclusters

fast_heatmap(data = results_overlap$gr_minus123_minusCumsum0.995$mclusters,
             data_type = "significant_uniq_data")

summary_overlap_results( results_overlap$gr_minus123_minusCumsum0.995$mclusters$significant_uniq_data)

# overlap analysis: clusters vs papers
analyze_gene_list_overlap(
  row_lists = filter_list_by_names( gr_database_blocked_gene_lists$marpiech_cluster_dex, c("cluster_18", "cluster_17"), operation = "!%in%"),
  col_lists =  set_data$gr_minus123_minusCumsum0.995$refine_gene_lists$gene_lists,
  reference_hgnc_vector = hgnc_symbols_vector_v110, 
  keep_original_data = TRUE, 
  fdr_threshold = 0.01, 
  overlap_threshold = 3
) -> results_overlap$gr_minus123_minusCumsum0.995$mclusters_minusUpDown

fast_heatmap(data = results_overlap$gr_minus123_minusCumsum0.995$mclusters_minusUpDown,
             data_type = "significant_uniq_data")

summary_overlap_results(results_overlap$gr_minus123_minusCumsum0.995$mclusters_minusUpDown$significant_uniq_data)


# overlap analysis: clusters vs papers
analyze_gene_list_overlap(
  row_lists = filter_list_by_names( gr_database_blocked_gene_lists$marpiech_cluster_dex, c("cluster_18", "cluster_17", "cluster_4"), operation = "!%in%"),
  col_lists =  set_data$gr_minus123_minusCumsum0.995$refine_gene_lists$gene_lists,
  reference_hgnc_vector = hgnc_symbols_vector_v110, 
  keep_original_data = TRUE, 
  fdr_threshold = 0.01, 
  overlap_threshold = 3
) -> results_overlap$gr_minus123_minusCumsum0.995$mclusters_minusUpDown4

fast_heatmap(data = results_overlap$gr_minus123_minusCumsum0.995$mclusters_minusUpDown4,
             data_type = "significant_uniq_data",overlap_threshold = 3)

summary_overlap_results(data = results_overlap$gr_minus123_minusCumsum0.995$mclusters_minusUpDown4$significant_uniq_data)

################################################################################
# minus123 and minusCumsum0.999

papers_data_preprocessing %>% 
  filter(!(treatment %in% treatment_to_remove)) %>% 
  filter(!(hgnc_symbol %in% set_data$gr_source$rare_gene_lists$freq_123)) %>% 
  filter(!(hgnc_symbol %in% set_data$gr_source$master_gene_lists$master_cumsum_0.999)) %>% 
  refine_gene_lists(data = .,
                    columns = c("source", "tissue", "cell", "treatment", "treatment_type", "regulation", "comparison", "environment"),
                    cumsum_thresholds = cumsum_thresholds,
                    freq_thresholds = freq_thresholds,
                    keep_column = c("system", "simple_tissue", "method")) -> set_data$gr_minus123_minusCumsum0.999


# overlap analysis: clusters vs papers
analyze_gene_list_overlap(
  row_lists = gr_database_blocked_gene_lists$marpiech_cluster_dex,
  col_lists =  set_data$gr_minus123_minusCumsum0.999$refine_gene_lists$gene_lists,
  reference_hgnc_vector = hgnc_symbols_vector_v110, 
  keep_original_data = TRUE, 
  fdr_threshold = 0.01, 
  overlap_threshold = 3
) -> results_overlap$gr_minus123_minusCumsum0.999$mclusters

fast_heatmap(data = results_overlap$gr_minus123_minusCumsum0.999$mclusters,
             data_type = "significant_uniq_data")

summary_overlap_results( results_overlap$gr_minus123_minusCumsum0.999$mclusters$significant_uniq_data)

# overlap analysis: clusters vs papers
analyze_gene_list_overlap(
  row_lists = filter_list_by_names( gr_database_blocked_gene_lists$marpiech_cluster_dex, c("cluster_18", "cluster_17"), operation = "!%in%"),
  col_lists =  set_data$gr_minus123_minusCumsum0.999$refine_gene_lists$gene_lists,
  reference_hgnc_vector = hgnc_symbols_vector_v110, 
  keep_original_data = TRUE, 
  fdr_threshold = 0.01, 
  overlap_threshold = 3
) -> results_overlap$gr_minus123_minusCumsum0.999$mclusters_minusUpDown

fast_heatmap(data = results_overlap$gr_minus123_minusCumsum0.999$mclusters_minusUpDown,
             data_type = "significant_uniq_data")

summary_overlap_results(results_overlap$gr_minus123_minusCumsum0.999$mclusters_minusUpDown$significant_uniq_data)


# overlap analysis: clusters vs papers
analyze_gene_list_overlap(
  row_lists = filter_list_by_names( gr_database_blocked_gene_lists$marpiech_cluster_dex, c("cluster_18", "cluster_17", "cluster_4"), operation = "!%in%"),
  col_lists =  set_data$gr_minus123_minusCumsum0.999$refine_gene_lists$gene_lists,
  reference_hgnc_vector = hgnc_symbols_vector_v110, 
  keep_original_data = TRUE, 
  fdr_threshold = 0.01, 
  overlap_threshold = 3
) -> results_overlap$gr_minus123_minusCumsum0.999$mclusters_minusUpDown4

fast_heatmap(data = results_overlap$gr_minus123_minusCumsum0.999$mclusters_minusUpDown4,
             data_type = "significant_uniq_data",overlap_threshold = 3)

summary_overlap_results(data = results_overlap$gr_minus123_minusCumsum0.999$mclusters_minusUpDown4$significant_uniq_data)
