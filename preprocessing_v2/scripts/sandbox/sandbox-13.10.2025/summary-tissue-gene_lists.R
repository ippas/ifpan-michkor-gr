
papers_data_preprocessing$simple_tissue %>% unique

papers_data_preprocessing %>% 
  filter(simple_tissue == "lung") %>% 
  .$time %>% unique

papers_data_preprocessing %>% 
  filter(simple_tissue == "lung") %>% summarize_gene_dataset_and_enrichr() -> tmp

tmp$summary_up$freq_genes_per_paper %>% 
  filter(freq >= 3) %>% 
  .$data %>% 
  unlist %>% 
  cat(sep = "\n")



tmp$summary_up$freq_genes_per_paper %>% filter(freq >= 3) %>% .$data %>% unlist

tmp$summary_down$freq_genes_per_paper %>% filter(freq >= 3) %>% .$data %>% unlist %>% cat(sep = "\n")


tmp$summary_down$freq_genes_per_paper

papers_data_preprocessing %>% 
  filter(simple_tissue == "kidney") %>% .$source %>% unique 


papers_data_preprocessing %>% 
  mutate(label = paste(label, regulation, sep = "_")) %>% 
  # dplyr::select(label, hgnc_symbol) %>% 
  group_by(label) %>%
  nest() %>% 
  mutate(n_genes = map(data, ~ .x$hgnc_symbol %>% unique %>% length)) %>% 
  unnest(n_genes) %>% 
  unnest(data) %>% 
  ungroup %>% .$n_genes %>% unique() %>% as.data.frame() %>% 
  set_colnames("Freq") %>% 
  filter(Freq >= 10) %>% 
  .$source %>% unique %>% 
  as.data.frame()  filter(simple_tissue == "blood") %>% .$label %>% unique()

papers_data_preprocessing %>% 
  filter(simple_tissue == "brain") %>% .$source %>% unique 



tmp_blood$summary_down$enrichr_summary_source$enrichr_3pub$CellMarker_2024 %>% head

papers_data_preprocessing %>%
  group_by(simple_tissue) %>%
  summarise(unique_sources = paste(unique(source), collapse = ", ")) %>%
  arrange(simple_tissue)


papers_data_preprocessing %>% 
  filter(simple_tissue == "blood") %>% summarize_gene_dataset_and_enrichr() -> tmp_blood

tmp_blood$summary_up$freq_genes_per_paper %>% 
  filter(freq >= 3) %>% .$data %>% unlist %>% 
  cat(sep = "\n")

  

papers_data_preprocessing %>% .$time %>% unique
