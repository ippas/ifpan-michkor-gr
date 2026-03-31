gene_list_log2ratio_top_50 %>% 
  select(-data) %>% 
  unnest(top_50_genes) %>% 
  .$hgnc_symbol %>% table %>% 
  as.data.frame() %>% 
  set_colnames(c("hgnc_symbol", "freq")) %>% 
  .$freq %>% 
  table %>% 
  as.data.frame() %>% 
  set_colnames(c("freq_genes", "n_genes")) %>% 
  mutate(propability = n_genes/number_genes) %>% 
  mutate(cumsum = cumsum(propability)) %>% 
  mutate(freq_genes = ifelse(cumsum > 0.97, "more_than_7", as.character(freq_genes))) %>% 
  group_by(freq_genes) %>% 
  summarise(n_genes = sum(n_genes), propability = sum(propability)) %>%
  ungroup() %>% 
  mutate(freq_genes = factor(freq_genes, levels = c(as.character(1:7), "more_than_7"))) %>% 
  arrange(freq_genes) %>% 
  mutate(cumsum = cumsum(propability)) -> data_pie_chart


################################################################################  
# dane do wykresu venna, dla wszystkich genów które mają wartości log2ratio
gene_list_log2ratio_top_50 %>% 
  select(-data) %>% 
  unnest(top_50_genes) %>% 
  ungroup %>% 
  select(hgnc_symbol, regulation) %>% 
  unique %>% 
  group_by(hgnc_symbol) %>% 
  nest() %>% 
  mutate(n_regulation = map(data, ~length(unique(.x$regulation)))) %>% 
  unnest() %>% 
  unique() %>% 
  mutate(group = ifelse(n_regulation == 2, "up_down", regulation)) %>% 
  select(hgnc_symbol, group) %>% 
  unique %>% 
  .$group %>% 
  table


