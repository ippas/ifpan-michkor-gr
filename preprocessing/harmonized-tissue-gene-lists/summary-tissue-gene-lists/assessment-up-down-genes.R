top50_genes_sum_rs %>% 
  select(hgnc_symbol, regulation) %>% 
  unique %>% 
  group_by(hgnc_symbol) %>% 
  nest() %>% 
  mutate(n_regulation = map(data, ~length(unique(.x$regulation)))) %>% 
  unnest() %>% 
  unique() %>% 
  mutate(group = ifelse(n_regulation == 2, "up_down", regulation)) %>% 
  select(hgnc_symbol, group) %>% 
  ungroup() %>% 
  unique() %>% 
  .$group %>% table

top50_genes_sum_rs %>% 
  select(hgnc_symbol, regulation) %>% 
  unique %>% 
  group_by(hgnc_symbol) %>% 
  nest() %>% 
  mutate(n_regulation = map(data, ~length(unique(.x$regulation)))) %>% 
  unnest() %>% 
  unique() %>% 
  mutate(group = ifelse(n_regulation == 2, "up_down", regulation)) %>% 
  select(hgnc_symbol, group) %>% 
  ungroup() %>% 
  filter(group == "up_down") %>% 
  unique %>% 
  .$hgnc_symbol -> double_reg_tissue_signatures

top50_genes_sum_rs %>% 
  filter(hgnc_symbol %in% double_reg_tissue_signatures) %>% 
  .$regulation %>% table


top50_genes_sum_rs %>% 
  filter(hgnc_symbol %in% double_reg_tissue_signatures) %>% 
  .$simple_tissue %>% 
  table


top50_genes_sum_rs %>% 
  filter(hgnc_symbol %in% double_reg_tissue_signatures) %>% 
  .$simple_tissue %>% 
  table

lapply(unique(top50_genes_sum_rs$simple_tissue), 
       function(x){
         top50_genes_sum_rs %>% 
           filter(hgnc_symbol %in% double_reg_tissue_signatures) %>% 
           filter(simple_tissue == x) %>% 
           group_by(hgnc_symbol) %>% 
           nest() %>% 
           mutate(n_reg = map(data, ~nrow(.x))) %>% 
           unnest(n_reg) %>% 
           filter(n_reg > 1) %>% 
           unnest(data) %>% unnest(data)
       })
       