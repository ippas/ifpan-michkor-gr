gr_genes_signatures_multi_approach_df %>% 
  filter(signature_name == "brain_up") %>% select(hgnc_symbol) %>% 
  print(index = F)
  


mdd2025_no23andMe_eur_gwas <- read.table("data/pgc-genes/pgc_mdd2025_no23andMe_eur.tsv", header = TRUE)


mdd2025_no23andMe_eur_gwas_preprocessing <- mdd2025_no23andMe_eur_gwas %>%
  mutate(
    gene_symbol = "None",
    signature_derivation = "depression_PMID:39814019",
    signature_name = "depression_PMID:39814019",
    range_plus_minus = "None"
  ) %>%
  select(
    rsID, chromosome, position, pvalue,
    gene_symbol, signature_name, signature_derivation, range_plus_minus
  )


protein_coding_intersection_data_GRCh37 <- read.table("data/pgc-genes/depression-PMID39914019/intersection_results_all_protein_coding_100000.tsv", header = TRUE)

################################################################################
protein_coding_intersection_data_GRCh37 %>% 
  filter(pvalue < 0.000001) %>% 
  .$gene_symbol %>% unique -> pmid39814019_1e06_100000


gr_genes_signatures_multi_approach_df %>% 
  filter(signature_derivation %in% c("universal", "tissue_cell")) %>% 
  select(-signature_derivation) %>% 
  group_by(signature_name) %>% 
  summarise(genes = list(hgnc_symbol), .groups = "drop") %>% 
  deframe() -> lists_overlap_depression

# large_gr_signatures$gr_full_signatures_by_tissue -> lists_overlap_depression


lists_overlap_depression$depression_pmid39814019_1e06_100000 <- pmid39814019_1e06_100000


# chi2_test_genesets(set1 = lists_overlap_depression$depression_pmid39814019_1e06_100000,
#                   set2 = lists_overlap_depression$`small-intestine_up`,
#                   alternative = "greater",
#                   background_genes = hgnc_symbols_vector_v110)
# 
# fisher_test_genes(gene_set_1 = lists_overlap_depression$depression_pmid39814019_1e06_100000,
#                   gene_set_2 = lists_overlap_depression$brain_glia_up,
#                    alternative = "less",
#                    background_genes = hgnc_symbols_vector_v110)


# protein_coding_intersection_data_GRCh37 %>% 
#   filter(pvalue < 0.000001) %>% 
#   filter(gene_symbol %in% {
#     gr_genes_signatures_multi_approach_df %>% 
#       filter(signature_name == "brain_up") %>% 
#       .$hgnc_symbol
#   }) %>% 
#   group_by(gene_symbol) %>% 
#   slice_min(pvalue, n = 1, with_ties = FALSE)


  

chi2_results_depression <- perform_chi2_tests(lists_overlap_depression, hgnc_symbols_vector_v110)

processing_overlap_results(data = chi2_results_depression ,
                           rows_to_filter = rownames(chi2_results_depression$p_value_matrix),
                           # rows_to_filter = tissues_clusters[10:27],
                           cols_to_filter = rownames(chi2_results_depression$p_value_matrix),
                           overlap_threshold = 0,
                           fdr_threshold = 1,
                           genes_list =lists_overlap_depression) -> depression_grSignatures_data

depression_grSignatures_data$original_data$df

depression_grSignatures_data$significant_data$df %>% 
  filter(Var2 == "depression_pmid39814019_1e06_100000") %>% 
  filter(gene_overlap_count > 0) %>% 
  arrange(p_value) %>% 
  select(-c(Var2, fdr))



fisher_test_genes <- function(gene_set_1, gene_set_2, background_genes = NULL, alternative = "greater") {
  # Jeśli nie podano tła, to jako unię wszystkich genów
  if (is.null(background_genes)) {
    background_genes <- union(gene_set_1, gene_set_2)
  }
  
  # Oblicz przecięcie
  overlap <- length(intersect(gene_set_1, gene_set_2))
  
  # Liczby do tabeli kontyngencji
  only_1 <- length(setdiff(gene_set_1, gene_set_2))
  only_2 <- length(setdiff(gene_set_2, gene_set_1))
  neither <- length(setdiff(background_genes, union(gene_set_1, gene_set_2)))
  
  # Tabela 2x2
  contingency_table <- matrix(c(overlap, only_1, only_2, neither), nrow = 2,
                              dimnames = list(GeneSet1 = c("Yes", "No"),
                                              GeneSet2 = c("Yes", "No")))
  
  # Test Fishera
  fisher_result <- fisher.test(contingency_table, alternative = alternative)
  
  return(list(
    p_value = fisher_result$p.value,
    odds_ratio = fisher_result$estimate,
    overlap_genes = intersect(gene_set_1, gene_set_2),
    contingency_table = contingency_table,
    result = fisher_result
  ))
}


chi2_test_genesets <- function(set1, set2, background_genes, alternative = "greater") {
  # Oblicz podstawowe wartości
  overlap <- length(intersect(set1, set2))
  only_set1 <- length(setdiff(set1, set2))
  only_set2 <- length(setdiff(set2, set1))
  
  # Liczba genów spoza obu list
  neither <- length(setdiff(background_genes, union(set1, set2)))
  
  # Macierz kontyngencji 2x2
  contingency_table <- matrix(
    c(overlap, only_set1,
      only_set2, neither),
    nrow = 2,
    dimnames = list(set1 = c("Yes", "No"),
                    set2 = c("Yes", "No"))
  )
  
  # Test chi-kwadrat
  chi2_result <- chisq.test(contingency_table, correct = FALSE)
  
  return(list(
    p_value = chi2_result$p.value,
    chi2_statistic = chi2_result$statistic,
    contingency_table = contingency_table,
    overlap_genes = intersect(set1, set2),
    alternative = alternative  # tylko informacyjnie
  ))
}

