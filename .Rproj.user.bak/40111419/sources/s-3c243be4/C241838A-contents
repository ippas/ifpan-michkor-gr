read.csv("data/phenotypes/panukbiobank-phenotype-category.csv") %>% 
  filter(category_name == "Depression") %>% 
  .$model_name %>% unique %>%  paste(., collapse = "|") -> depression_models_pattern


# List all files in the directory with full names
list.files("data/prs-models-pan-biobank-uk/", full.names = TRUE) -> model_files_vector

model_files_vector %>% 
  keep(~ grepl(depression_models_pattern, .x)) %>% 
  set_names(map(., ~ basename(.x) %>% tools::file_path_sans_ext())) %>%
  # Read YAML content from each file
  map(read_yaml) %>%
  # Extract genes information and name the list elements after the files
  map(~ .x$description$genes) %>% 
  lapply(., unique) -> depression_category_genes_list

######################################################################Sys.time() -> start_time
permutation_mental_health_results2 <- perform_overlap_permutation_analysis_parallel(
  permutations = 100,  # Number of permutations
  seed = 123,          # Initial seed for reproducibility
  num_cores = 35,       # Number of cores to use
  # Additional arguments for analyze_random_gene_sets
  reference_hgnc_vector = hgnc_symbols_vector_v110,
  size_reference_list = papers_gene_list[!grepl("marpiech", names(papers_gene_list))],
  comparison_gene_list = depression_category_genes_list,
  overlap_threshold = 3,
  fdr_threshold = 0.01
)

Sys.time() -> end_time

start_time - end_time##########
# analysis
phenotypes_biobank_genes_list <- icd10f_genes_list

papers_gene_list %>% names

# papers_gene_list[!grepl("marpiech", names(papers_gene_list))]

c(papers_gene_list, phenotypes_biobank_genes_list) %>% 
  lapply(., unique) -> phenotypes_clusters_list 

chi2_results_phenotypes <- perform_chi2_tests(phenotypes_clusters_list, hgnc_symbols_vector_v110)

processing_overlap_results(data = chi2_results_phenotypes,
                           rows_to_filter = !rownames(chi2_results_phenotypes$p_value_matrix) %in% names(papers_gene_list),
                           # rows_to_filter = tissues_clusters[10:27],
                           cols_to_filter = names(papers_gene_list),
                           overlap_threshold = 3,
                           fdr_threshold = 0.05,
                           genes_list = phenotypes_clusters_list) -> clusters_phenotypes_data


# clusters_phenotypes_data$significant_uniq_data$df %>% 
#   mutate(Var1 = str_replace_all(Var1, "biobankuk.*both_sexes-", "")) %>% 
#   mutate(Var1 = str_replace_all(Var1, "-EUR-1e-05", "")) %>% select(-number_overlap) %>% 
#   arrange(fdr) %>% 
#   write.csv(., file = "results/tables/psychiatric-papers-overlap-fdr0.01.csv", quote = F, row.names = F)


clusters_phenotypes_data$significant_uniq_data$rows %>% gsub("biobankuk.*both_sexes-", "", .) -> phenotyepes_mapping

names(phenotyepes_mapping) <- clusters_phenotypes_data$significant_uniq_data$rows 


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
  # col_mapping_vector =  clusters_mapping,
  row_mapping_vector =  phenotyepes_mapping,
  fdr_threshold = 0.01,
  fdr_thresholds = c(0.05, 0.0001),
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
  row_names_gp = gpar(fontsize = 12),
  column_names_gp = gpar(fontsize = 12),
  column_names_rot = 45,
  column_names_side = "top",
  overlap_threshold = 3
)


clusters_phenotypes_data$significant_uniq_data$df %>% arrange(fdr) %>% head(8) %>% 
  .$overlap_genes %>% as.list() %>% 
  lapply(., function(x){unlist(strsplit(x, ","))}) %>% 
  unlist()

clusters_phenotypes_data$significant_uniq_data$df %>% arrange(fdr) %>% 
  write.table(., "results/tables/psychiatric-overlap/mental-health-overlap.tsv", row.names = F, col.names = T, quote = F, sep = "\t")
