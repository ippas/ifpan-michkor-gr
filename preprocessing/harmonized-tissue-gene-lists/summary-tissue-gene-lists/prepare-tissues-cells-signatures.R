# tissue_cells

rbind(top50_per_tissue, top50_per_celltype) %>%
  unnest(data) %>%
  group_by(simple_tissue, regulation, hgnc_symbol) %>%
  summarise(sum_rs = sum(rank_score, na.rm = TRUE), .groups = "drop") %>%
  group_by(simple_tissue, regulation) %>%
  arrange(desc(sum_rs)) %>%
  slice_head(n = 50) %>% 
  mutate(signature_name = paste(simple_tissue, regulation, sep = "_")) 

  
  

# universal
top50_genes_sum_rs %>% 
  unnest(data) %>%
  group_by(simple_tissue, regulation, hgnc_symbol) %>%
  summarise(sum_rs = sum(rank_score, na.rm = TRUE), .groups = "drop") %>%
  group_by(simple_tissue, regulation) %>%
  arrange(desc(sum_rs)) %>%
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
  slice_max(second_sum_rs, n = 50) %>% 
  select(regulation, hgnc_symbol) %>% 
  ungroup %>% 
  mutate(signature_name = paste0("universal_", regulation)) %>% 
  select(hgnc_symbol, signature_name)
  

top50_genes_sum_rs %>% 
  unnest(data) %>%
  group_by(simple_tissue, regulation, hgnc_symbol) %>%
  summarise(sum_rs = sum(rank_score, na.rm = TRUE), .groups = "drop") %>%
  group_by(simple_tissue, regulation) %>%
  arrange(desc(sum_rs)) %>%
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
  slice_max(second_sum_rs, n = 50) %>% 
  select(-data) %>% 
  filter(regulation == "up") %>% .$second_sum_rs %>% mean
