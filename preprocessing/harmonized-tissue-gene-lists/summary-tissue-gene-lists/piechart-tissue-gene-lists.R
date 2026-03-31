top50_genes_sum_rs %>% 
  select(-data) %>% 
  .$hgnc_symbol %>% 
  table %>% 
  as.data.frame() %>% 
  set_colnames(c("hgnc_symbol", "freq")) %>% 
  .$freq %>% 
  table %>% 
  as.data.frame() %>% 
  set_colnames(c("freq_genes", "n_genes")) %>% 
  mutate(propability = n_genes/1055) %>% 
  mutate(cumsum = cumsum(propability)) %>% 
  mutate(freq_genes = ifelse(cumsum > 0.995, "more_than_4", as.character(freq_genes))) %>% 
  group_by(freq_genes) %>% 
  summarise(n_genes = sum(n_genes), propability = sum(propability)) %>%
  ungroup() %>% 
  mutate(freq_genes = factor(freq_genes, levels = c(as.character(1:7), "more_than_4"))) %>% 
  arrange(freq_genes) %>% 
  mutate(cumsum = cumsum(propability)) -> data_pie_chart

svg("results/figures/harmonized-tissue-gene-lists/summary-tissue-gene-lists/piechart-tissue-gene-lists.svg", width = 8, height = 8)
ggplot(data_pie_chart, aes(x = "", y = n_genes, fill = freq_genes)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_brewer(palette = "Pastel2") +  # Użycie poprawnej funkcji dla palety RColorBrewer
  labs(fill = "Gene Frequency") +
  ggtitle("Distribution of Gene Frequency")
dev.off()

# piechart for down genes

top50_genes_sum_rs %>% 
  filter(regulation == "down") %>% 
  select(-data) %>% 
  .$hgnc_symbol %>% 
  table %>% 
  as.data.frame() %>% 
  set_colnames(c("hgnc_symbol", "freq")) %>% 
  .$freq %>% 
  table %>% 
  as.data.frame() %>% 
  set_colnames(c("freq_genes", "n_genes")) %>% 
  mutate(propability = n_genes/625) %>% 
  mutate(cumsum = cumsum(propability)) ->  data_pie_chart


svg("results/figures/harmonized-tissue-gene-lists/summary-tissue-gene-lists/piechart-tissue-gene-lists-down.svg", width = 8, height = 8)
ggplot(data_pie_chart, aes(x = "", y = n_genes, fill = freq_genes)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_brewer(palette = "Pastel2") +  # Użycie poprawnej funkcji dla palety RColorBrewer
  labs(fill = "Gene Frequency") +
  ggtitle("Distribution of Gene Frequency")
dev.off()

top50_genes_sum_rs %>% 
  filter(regulation == "up") %>% 
  select(-data) %>% 
  .$hgnc_symbol %>% 
  table %>% 
  as.data.frame() %>% 
  set_colnames(c("hgnc_symbol", "freq")) %>% 
  .$freq %>% 
  table %>% 
  as.data.frame() %>% 
  set_colnames(c("freq_genes", "n_genes")) %>% 
  mutate(propability = n_genes/487) %>% 
  mutate(cumsum = cumsum(propability)) %>% 
  mutate(freq_genes = ifelse(cumsum > 0.99, "more_than_4", as.character(freq_genes))) %>% 
  group_by(freq_genes) %>% 
  summarise(n_genes = sum(n_genes), propability = sum(propability)) %>%
  ungroup() %>% 
  mutate(freq_genes = factor(freq_genes, levels = c(as.character(1:7), "more_than_4"))) %>% 
  arrange(freq_genes) %>% 
  mutate(cumsum = cumsum(propability)) -> data_pie_chart

svg("results/figures/harmonized-tissue-gene-lists/summary-tissue-gene-lists/piechart-tissue-gene-lists-up.svg", width = 8, height = 8)
ggplot(data_pie_chart, aes(x = "", y = n_genes, fill = freq_genes)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_brewer(palette = "Pastel2") +  # Użycie poprawnej funkcji dla palety RColorBrewer
  labs(fill = "Gene Frequency") +
  ggtitle("Distribution of Gene Frequency")
dev.off()
