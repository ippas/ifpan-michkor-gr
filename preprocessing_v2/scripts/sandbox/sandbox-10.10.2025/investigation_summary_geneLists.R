AllBrainGeneDf %>%
  filter(n_genes >= 10) %>%
  filter(!(time %in% c("240h", "720h", "240h_vs_720h", "672h", "3months", "3weeks", "10days"))) %>%
  group_by(label) %>%
  nest() %>%
  mutate(n_clusterK = map_int(data, ~ sum(.x$hgnc_symbol %in% AllBrain2BrainSignatures2GlobalMarpiechClusters$cluster_K))) %>%
  filter(n_clusterK >= 5) %>%
  unnest(data) -> tmp_dataset

dataset_list$n10_clusterK_min5_short_time %>% .$source %>% unique

tmp <- summarize_gene_dataset(df = AllBrainGeneDf, dataset_name = "test")
tmp



tmp <- summarize_gene_dataset(df = tmp_dataset, dataset_name = "test")
tmp


tmp_dataset %>% 
  select(label, hgnc_symbol) %>% 
  unique %>% .$hgnc_symbol %>% table %>% as.data.frame() %>% .$Freq %>% table


AllBrainGeneDf %>% 
  select(label, hgnc_symbol) %>% 
  unique() %>% 
  group_by(label) %>% 
  nest() %>% 
  mutate(n_genes = map(data, ~ .x %>% nrow)) %>% 
  unnest(n_genes) %>% .$n_genes %>% summary
  


AllBrainGeneDf %>% 
  filter(n_genes)
