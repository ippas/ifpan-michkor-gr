
data_filtered %>% 
  filter(signature_name == "michkor_indicate_T.S1_Cx43DGE_n67") %>% 
  select(rsID, pvalue, gene_symbol) %>% 
  group_by(gene_symbol) %>% 
  mutate(n_rsID = n()) %>% 
  slice_min(pvalue, with_ties = FALSE)
  # gene-overlap-slezak-vs-depressionGWAS
  
slezak_gene_signatures_preprocessing

gene_lists <- slezak_signatures_list$michkor_indicate_T.S1_Cx43DGE

gene_lists <- list()

for (i in 1:1000) {
  random_genes <- sample(protein_coding_intersection_data_GRCh37$gene_symbol, 62, replace = FALSE)
  gene_lists[[paste0("random", i)]] <- random_genes
}

gene_lists <- list()
pb <- progress::progress_bar$new(
  format = "  Generating gene lists [:bar] :percent eta: :eta",
  total = 1000, clear = FALSE, width = 60
)

for (i in 1:1000) {
  random_genes <- sample(protein_coding_intersection_data_GRCh37$gene_symbol, 67, replace = FALSE)
  gene_lists[[paste0("random", i)]] <- random_genes
  pb$tick()
}


# gene_lists <- gr_genes_signatures_multi_approach_list

protein_coding_intersection_data_GRCh37 %>% 
  filter(pvalue < 1e-4) %>% 
  .$gene_symbol %>% unique %>% length()


depression_GWAS_PMID39814019_1e4 <- 
  protein_coding_intersection_data_GRCh37 %>% 
  filter(pvalue < 1e-4) %>%
  .$gene_symbol %>% unique 

# gene_lists$michkor_indicate_T.S1_Cx43DGE_n62 <-


chi2_results_slezak <- perform_chi2_tests(
  c(slezak_signatures_list$michkor_indicate_T.S1_Cx43DGE, gene_lists),
  {protein_coding_intersection_data_GRCh37$gene_symbol %>% unique},
  verbose = F)


chi2_results_slezak$number_overlap_matrix %>% tail

{chi2_results_slezak$number_overlap_matrix["depression_GWAS_PMID39814019_1e4",] > 15} %>% sum

chi2_results_slezak$number_overlap_matrix %>% colnames() %>% [1:10]

processing_overlap_results(data = chi2_results_slezak,
                           rows_to_filter = names(gene_lists),
                           cols_to_filter = c("michkor_indicate_T.S1_Cx43DGE_n67"),
                           overlap_threshold = 0,
                           fdr_threshold = 1,
                           genes_list = gene_lists) -> slezak_depression_GWAS_depression_PMID39814019_preprocessing

slezak_depression_GWAS_depression_PMID39814019_preprocessing$significant_uniq_data$df %>% head
  filter(Var2 == "michkor_indicate_T.S1_Cx43DGE_n67") %>% select(-fdr) %>% head 
  # filter(Var1 != "michkor_indicate_T.S1_Cx43DGE_n67") %>% 
  filter(grepl("random", Var1)) %>% 
  .$gene_overlap_count %>% hist


chi2_results_slezak$p_value_matrix["depression_GWAS_PMID39814019_1e4", ]
chi2_results_slezak$overlap_genes_matrix


protein_coding_intersection_data_GRCh37[sample(nrow(protein_coding_intersection_data_GRCh37), 67), ] %>% 
  .$gene_symbol

protein_coding_intersection_data_GRCh37 %>% 
  sample(67)
