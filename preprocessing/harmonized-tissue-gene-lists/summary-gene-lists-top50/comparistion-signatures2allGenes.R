top50_genes_sum_rs %>% 
  select(-data) %>% 
  .$hgnc_symbol %>% table %>% 
  as.data.frame() %>% 
  set_colnames(c("hgnc_symbol", "freq")) %>% 
  arrange(desc(freq)) -> freq_top_genes

papers_data_preprocessing %>%
  .$hgnc_symbol %>%
  table() %>%
  as.data.frame() %>% 
  set_colnames(c("hgnc_symbol", "freq")) -> freq_all_genes

number_genes <- freq_top_genes$hgnc_symbol %>% unique() %>% length()

freq_all_genes %>% 
  filter(hgnc_symbol %in% freq_top_genes$hgnc_symbol) %>%
  count(freq) %>%  # liczy ile genów ma daną wartość freq
  set_colnames(c("freq_genes", "n_genes")) %>% 
  mutate(
    propability = n_genes / number_genes,
    cumsum = cumsum(propability),
    rev_cumsum = 1 - cumsum,
    freq_genes = ifelse(freq_genes > 12, "more_than_12", as.character(freq_genes))
  ) %>% 
  group_by(freq_genes) %>% 
  summarise(
    n_genes = sum(n_genes),
    propability = sum(propability),
    .groups = "drop"
  ) %>%
  mutate(
    freq_genes = factor(freq_genes, levels = c(as.character(1:12), "more_than_12"))
  ) %>% 
  arrange(freq_genes) %>% 
  mutate(cumsum = cumsum(propability))
