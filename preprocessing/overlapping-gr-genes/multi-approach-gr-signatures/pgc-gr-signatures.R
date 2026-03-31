

# gwas pgc
c(categorized_gene_lists$phenotypes_PGC$gene_lists, 
  gr_genes_signatures_multi_approach_list) %>% 
  lapply(., unique) -> gene_lists


chi2_results_phenotypesPGC_signatures <- perform_chi2_tests(
  gene_lists, 
  hgnc_symbols_vector_v110)

###
processing_overlap_results(data = chi2_results_phenotypesPGC_signatures,
                           rows_to_filter = names(gr_genes_signatures_multi_approach_list),
                           cols_to_filter = names(categorized_gene_lists$phenotypes_PGC$gene_lists),
                           overlap_threshold = 3,
                           fdr_threshold = 1,
                           genes_list = gene_lists) -> pgc_signatures_gr


pgc_signatures_gr$significant_data$df %>% 
  filter(p_value < 0.05) %>% 
  rename(Var1 = "signature_name") %>% 
  rename(Var2 = "phenotype") %>% 
  write.table("results/gr-signatures/overlap-PGC-tissuCellUniversalJuszczak-23.03.2025.tsv",
              quote = FALSE,
              row.names = FALSE,
              col.names = TRUE,
              sep = "\t")

# combine diseases pgc
c(categorized_gene_lists$diseases_PGC$gene_lists, 
  gr_genes_signatures_multi_approach_list) %>% 
  lapply(., unique) -> gene_lists


chi2_results_diseasesPGC_signatures <- perform_chi2_tests(
  gene_lists, 
  hgnc_symbols_vector_v110)

###
processing_overlap_results(data = chi2_results_diseasesPGC_signatures,
                           cols_to_filter = names(gr_genes_signatures_multi_approach_list),
                           rows_to_filter = names(categorized_gene_lists$diseases_PGC$gene_lists),
                           overlap_threshold = 3,
                           fdr_threshold = 1,
                           genes_list = gene_lists) -> diseases_pgc_signatures_gr


diseases_pgc_signatures_gr$significant_uniq_data$df %>% 
  filter(p_value < 0.05)

pgc_signatures_gr$significant_uniq_data$df$overlap_genes %>% 
  lapply(., str_split, ",") %>%  unlist %>% 
  unique() %>% length()
