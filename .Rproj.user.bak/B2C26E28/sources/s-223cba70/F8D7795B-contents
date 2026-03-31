top50_genes_sum_rs %>% 
  nest %>% 
  mutate(data = map(data, ~ .x %>% 
                       arrange(desc(sum_rs)) %>% 
                       mutate(second_rs = row_number(sum_rs)))) %>%
  unnest(data) %>% 
  ungroup() %>% 
  group_by(hgnc_symbol, regulation) %>% 
  nest() %>% 
  mutate(second_sum_rs = map(data, ~ .x %>% .$second_rs %>% sum)) %>% 
  unnest(second_sum_rs) %>% 
  arrange(desc(second_sum_rs)) %>% 
  ungroup() %>% 
  group_by(regulation) %>% 
  slice_max(second_sum_rs, n = 50) -> universal_gr_genes_nest

universal_gr_genes_nest %>% 
  ungroup %>% 
  # mutate(list_name = ifelse(regulation == "up", "up", "down")) %>% 
  unnest(data) %>% 
  unnest(data) %>% 
  ungroup -> universal_gr_genes
  
universal_gr_genes %>% filter(hgnc_symbol == "ZBTB16")  %>% 
  select(regulation, simple_tissue) %>% 
  unique()
  
universal_gr_genes$hgnc_symbol %>% unique() %>% length()


universal_gr_genes_nest %>% filter(regulation == "down") %>%
  select(c(hgnc_symbol, second_sum_rs)) %>% 
  ungroup %>% 
  select(c(hgnc_symbol, second_sum_rs)) %>% as.data.frame()


# summary lists
universal_gr_genes_nest %>% 
  filter(regulation == "down") %>% 
  unnest(data) %>% 
  ungroup() %>% 
  group_by(hgnc_symbol) %>% 
  nest %>% 
  mutate(freq = map(data, ~nrow(.x))) %>% 
  unnest(freq) %>% 
  arrange(freq) %>% 
  filter(freq == 1) %>% unnest(data)
  .$hgnc_symbol %>% 
  cat(sep = "\n")


universal_gr_genes_nest$hgnc_symbol 


papers_data_preprocessing %>% 
  .$hgnc_symbol %>% 
  table %>% 
  as.data.frame() %>% 
  set_colnames(c("hgnc_symbol", "freq")) %>% 
  filter(hgnc_symbol %in% universal_gr_genes_nest$hgnc_symbol ) %>% 
  arrange(freq) %>% 
  filter(freq >  12)
