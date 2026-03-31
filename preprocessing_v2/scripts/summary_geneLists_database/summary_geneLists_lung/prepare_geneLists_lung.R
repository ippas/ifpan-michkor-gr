papers_data_preprocessing %>% 
  mutate(label = paste(label, regulation, sep = "_")) %>% 
  group_by(label) %>%
  nest() %>% 
  mutate(n_genes = map(data, ~ .x$hgnc_symbol %>% unique %>% length)) %>% 
  unnest(n_genes) %>% 
  unnest(data) %>% 
  ungroup %>% 
  as.data.frame() -> papers_data_preprocessing


papers_data_preprocessing %>% 
  filter(simple_tissue == "lung") -> AllLungGeneDf
