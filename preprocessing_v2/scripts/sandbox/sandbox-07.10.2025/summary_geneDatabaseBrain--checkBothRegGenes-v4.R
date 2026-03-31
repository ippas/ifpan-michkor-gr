# ##############################################################################
# ---- uses data ----
# ##############################################################################
summary_table 
gene_summary_list

dataset_list$n10_short_time %>% 
  filter(hgnc_symbol == "MAP3K6")



gene_summary_list$all$regulation_summary %>% 
  filter(regulation_type == "only_up") %>% 
  .$genes %>% .[[1]] %>% 
  as.data.frame() %>% 
  set_colnames(c("hgnc_symbol")) %>% 
  filter(hgnc_symbol %in% AllBrain2BrainSignatures2GlobalMarpiechClusters$cluster_P)


summary_table %>% 
  as.data.frame() %>% 
  mutate(percent_both = both_n / n_genes*100) %>% 
  mutate(percent_unique = n_genes_freq1 / n_genes*100)

