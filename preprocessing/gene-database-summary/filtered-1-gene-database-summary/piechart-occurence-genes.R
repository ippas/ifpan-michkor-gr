# number records
filt_1_gr_database %>%  nrow 

# number genes
filt_1_gr_database %>% .$hgnc_symbol %>% unique() %>% length() -> number_genes


filt_1_gr_database %>%   
  select(hgnc_symbol, hgnc_occurence) %>%
  unique() %>%
  arrange(desc(hgnc_occurence)) %>% 
  .$hgnc_occurence %>% 
  table %>% 
  as.data.frame() %>% 
  set_colnames(c("freq_genes", "n_genes")) %>% 
  mutate(propability = n_genes/number_genes) %>% 
  mutate(cumsum = cumsum(propability)) %>% 
  mutate(freq_genes = ifelse(cumsum > 0.94, "more_than_12", as.character(freq_genes))) %>%
  group_by(freq_genes) %>% 
  summarise(n_genes = sum(n_genes), propability = sum(propability)) %>%
  ungroup() %>%
  mutate(freq_genes = factor(freq_genes, levels = c(as.character(1:12), "more_than_12"))) %>%
  arrange(freq_genes) %>% 
  mutate(cumsum = cumsum(propability)) -> data_pie_chart


n <- 12
custom_palette <-c(brewer.pal(n = n + 1, name = "Paired")[-1], "white")

svg("results/figures/gene-database-summary/filtered-1-gene-database-summary/piechart-filt-1-genes-occurence.svg", width = 8, height = 8)
ggplot(data_pie_chart, aes(x = "", y = n_genes, fill = freq_genes)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_manual(values = custom_palette) +
  labs(fill = "Gene Frequency") +
  ggtitle("Distribution of Gene Frequency")
dev.off()
  
  
