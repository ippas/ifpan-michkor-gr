c(categorized_gene_lists$metabolome_gene_lists$gene_lists, 
  gr_genes_signatures_multi_approach_list) %>% 
  lapply(., unique) -> gene_lists


chi2_results_metabolone_signatures <- perform_chi2_tests(
  gene_lists, 
  hgnc_symbols_vector_v110)

###
chi2_results_phenotypes_signatures$number_overlap_matrix

###
processing_overlap_results(data = chi2_results_metabolone_signatures,
                           rows_to_filter = names(gr_genes_signatures_multi_approach_list),
                           cols_to_filter = names(categorized_gene_lists$metabolome_gene_lists$gene_lists),
                           overlap_threshold = 2,
                           fdr_threshold = 0.1,
                           genes_list = gene_lists) -> tmp


tmp$significant_data$df %>% 
  arrange(Var1) %>% 
  unique()

tmp$gene_list_sizes


################################################################################
c(categorized_gene_lists$nightingale_gene_lists$gene_lists, 
  gr_genes_signatures_multi_approach_list) %>% 
  lapply(., unique) -> gene_lists


chi2_results_metabolone_signatures <- perform_chi2_tests(
  gene_lists, 
  hgnc_symbols_vector_v110)

###
chi2_results_phenotypes_signatures$number_overlap_matrix

###
processing_overlap_results(data = chi2_results_metabolone_signatures,
                           rows_to_filter = names(gr_genes_signatures_multi_approach_list),
                           cols_to_filter = names(categorized_gene_lists$nightingale_gene_lists$gene_lists),
                           overlap_threshold = 3,
                           fdr_threshold = 0.01,
                           genes_list = gene_lists) -> tmp


tmp$significant_data$df %>% 
  arrange(Var1) %>% 
  unique() %>% 
  filter(Var1 == "kidney_up") %>%
  .$Var2 %>% unique()


# można wziąć, geny z tkanki i sprawdzić czy są regulowane przez GR w innych tkankach.