gene_list_log2ratio_min_50_genes %>% 
  str


gene_list_log2ratio_min_50_genes %>% 
  group_by(label_regulation) %>% 
  nest() %>% 
  mutate(top_50_genes = map(data, ~ .x %>% 
                              arrange(desc(abs(log2ratio))) %>% 
                              head(50) %>% 
                              mutate(rank_score = row_number(abs(log2ratio))))) -> gene_list_log2ratio_top_50

number_genes <- gene_list_log2ratio_top_50 %>% 
  select(-data) %>% 
  unnest(top_50_genes) %>% 
  .$hgnc_symbol %>% 
  unique() %>% 
  length()

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
  


svg("results/figures/harmonized-tissue-gene-lists/piechart-gene-list-top50.svg", width = 8, height = 8)
ggplot(data_pie_chart, aes(x = "", y = n_genes, fill = freq_genes)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_brewer(palette = "Pastel2") +  # Użycie poprawnej funkcji dla palety RColorBrewer
  labs(fill = "Gene Frequency") +
  ggtitle("Distribution of Gene Frequency")

dev.off()

rm(data_pie_chart)


gene_list_log2ratio_top_50 %>% 
  select(-data) %>% 
  
  
  filter(simple_tissue %in% c("lung", "blood", "brain"))


