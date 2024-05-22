papers_data_preprocessing  %>%  filter(source != "marpiech_tissues") %>% extract_keys_values(., "info", keys = "method") -> papers_data_preprocessing
columns <- c("source")
cumsum_thresholds <- c(0.99, 0.995, 0.996, 0.997, 0.998, 0.999)
freq_thresholds <- list(1, 2, 3, c(1, 2), c(1, 2, 3))

treatment_to_remove <- c("vehicle-DMSO", "vitamin-d3", "formoterol", "LPS", "TNFalpha", "PHA", "TNF", "vehicle-ethanol")

results_overlap <- list()

set_data <- list()

refine_gene_lists(data = papers_data_preprocessing,
                  columns = columns,
                  cumsum_thresholds = cumsum_thresholds,
                  freq_thresholds = freq_thresholds) -> set_data$gr_source_all

papers_data_preprocessing %>% 
  filter(!(treatment %in% treatment_to_remove)) %>% 
  refine_gene_lists(data = .,
                  columns = columns,
                  cumsum_thresholds = cumsum_thresholds,
                  freq_thresholds = freq_thresholds) -> set_data$gr_source



refine_gene_lists(data = filter(papers_data_preprocessing, regulation == "up"),
                  columns = c("source", "regulation"),
                  cumsum_thresholds = cumsum_thresholds,
                  freq_thresholds = freq_thresholds) -> set_data$gr_source_up

refine_gene_lists(data = filter(papers_data_preprocessing, regulation == "down"),
                  columns = c("source", "regulation"),
                  cumsum_thresholds = cumsum_thresholds,
                  freq_thresholds = freq_thresholds) -> set_data$gr_source_down

papers_data_preprocessing %>%
  filter(!(treatment %in% treatment_to_remove)) %>%
  refine_gene_lists(data = .,
                  columns = c("source", "tissue", "cell", "treatment", "treatment_type", "regulation", "comparison", "environment"),
                  cumsum_thresholds = cumsum_thresholds,
                  freq_thresholds = freq_thresholds,
                  keep_column = c("system", "simple_tissue", "method")) -> set_data$gr_list

################################################################################
set_data$gr_source$genes_distribution$genes_frequency %>% 
  .[, c(1:4)] %>% as.data.frame()



################################################################################
# simple informations
set_data$gr_list$refine_gene_lists$metadata$environment %>% table

set_data$gr_list$refine_gene_lists$metadata$treatment_type %>% table

set_data$gr_list$refine_gene_lists$metadata$treatment %>% unique() %>% cat(sep = "\n")

################################################################################
# investigation of size gene lists


set_data$gr_list$refine_gene_lists$metadata %>% dim

set_data$gr_list$refine_gene_lists$metadata %>% 
  ggplot(aes(x = gene_count)) + 
  geom_histogram() +
  theme_minimal()

set_data$gr_list$refine_gene_lists$metadata %>% 
  ggplot(aes(x = gene_count)) + 
  geom_histogram() +
  facet_wrap(~ method) +
  theme_minimal()

set_data$gr_list$refine_gene_lists$metadata %>% 
  filter(regulation %in% c("up", "down")) %>% 
  ggplot(aes(x = gene_count)) + 
  geom_histogram() +
  facet_wrap(~ regulation) +
  theme_minimal()

set_data$gr_list$refine_gene_lists$metadata %>% 
  mutate(decile = ntile(gene_count, 10))
    
set_data$gr_list$refine_gene_lists$metadata %>%
  mutate(decile = ntile(gene_count, 10)) %>%
  group_by(decile) %>%
  summarize(
    mean_gene_count = mean(gene_count),
    median_gene_count = median(gene_count),
    min_gene_count = min(gene_count),
    max_gene_count = max(gene_count),
    sd_gene_count = sd(gene_count)
  )  


set_data$gr_list$refine_gene_lists$metadata %>% 
  filter(is.na(method))

set_data$gr_list$refine_gene_lists$metadata %>% 
  # filter(source == "pmid:23110767")
  filter(tissue == "kidney")



################################################################################
# overlap analysis: clusters vs papers
analyze_gene_list_overlap(
  row_lists = gr_database_blocked_gene_lists$marpiech_cluster_dex,
  col_lists = set_data$gr_list$refine_gene_lists$gene_lists,
  reference_hgnc_vector = hgnc_symbols_vector_v110, 
  keep_original_data = TRUE, 
  fdr_threshold = 0.01, 
  overlap_threshold = 3
) -> results_overlap$raw_list$mcluster

fast_heatmap(data =   results_overlap$raw_list$mcluster,
             data_type = "significant_uniq_data")

summary_overlap_results(results_overlap$raw_list$mcluster$significant_uniq_data)


# overlap analysis: clusters (minus UP and DOWN) vs papers 
analyze_gene_list_overlap(
  row_lists = filter_list_by_names( gr_database_blocked_gene_lists$marpiech_cluster_dex, c("cluster_18", "cluster_17"), operation = "!%in%"),
  col_lists = set_data$gr_list$refine_gene_lists$gene_lists,
  reference_hgnc_vector = hgnc_symbols_vector_v110, 
  keep_original_data = TRUE, 
  fdr_threshold = 0.01, 
  overlap_threshold = 3
) -> results_overlap$raw_list$mcluster_minusUpDown

fast_heatmap(data =   results_overlap$raw_list$mcluster_minusUpDown,
             data_type = "significant_uniq_data")

summary_overlap_results(results_overlap$raw_list$mcluster_minusUpDown$significant_uniq_data)

# overlap analysis: clusters (minus UP and DOWN) vs papers 
analyze_gene_list_overlap(
  row_lists = filter_list_by_names( gr_database_blocked_gene_lists$marpiech_cluster_dex, c("cluster_18", "cluster_17"), operation = "!%in%"),
  col_lists = set_data$gr_list$refine_gene_lists$gene_lists,
  reference_hgnc_vector = hgnc_symbols_vector_v110, 
  keep_original_data = TRUE, 
  fdr_threshold = 0.01, 
  overlap_threshold = 10
) -> results_overlap$raw_list$mcluster_minusUpDown_overlap10

fast_heatmap(data =   results_overlap$raw_list$mcluster_minusUpDown_overlap10,
             data_type = "significant_uniq_data",
             overlap_threshold = 10)

summary_overlap_results(results_overlap$raw_list$mcluster_minusUpDown_overlap10$significant_uniq_data)

restore_overlap_analysis$mcluster_vs_full_lists$significant_uniq_data$df$gene_overlap_count %>% 
  summary 

restore_overlap_analysis$mcluster_vs_full_lists$significant_uniq_data$df %>% 
  filter(!(Var2 %in% c("cluster_17", "cluster_18"))) %>% 
  .$gene_overlap_count %>% 
  summary 

restore_overlap_analysis$mcluster_vs_full_lists$significant_uniq_data$df %>% 
  filter(Var2 %in% c("cluster_17", "cluster_18")) %>% 
  .$gene_overlap_count %>% 
  summary 


restore_overlap_analysis$mcluster_vs_full_lists$significant_uniq_data$df %>% 
  filter(!(Var2 %in% c("cluster_17", "cluster_18"))) %>% 
  filter(gene_overlap_count > 8) %>% 
  filter(Var2 == "cluster_4")
  


