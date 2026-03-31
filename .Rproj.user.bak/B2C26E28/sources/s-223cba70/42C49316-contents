# number records
filt_1_12_gr_database %>%  nrow 

# number genes
filt_1_12_gr_database %>% .$hgnc_symbol %>% unique() %>% length() -> number_genes


filt_1_12_gr_database %>%   
  select(hgnc_symbol, hgnc_occurence) %>%
  unique() %>%
  arrange(desc(hgnc_occurence)) %>% 
  .$hgnc_occurence %>% 
  table %>% 
  as.data.frame() %>% 
  set_colnames(c("freq_genes", "n_genes")) %>% 
  mutate(propability = n_genes/number_genes) %>% 
  mutate(cumsum = cumsum(propability)) %>% 
  mutate(freq_genes = ifelse(cumsum > 0.84, "more_than_22", as.character(freq_genes))) %>%
  group_by(freq_genes) %>% 
  summarise(n_genes = sum(n_genes), propability = sum(propability)) %>%
  ungroup() %>%
  mutate(freq_genes = factor(freq_genes, levels = c(as.character(13:22), "more_than_22"))) %>%
  arrange(freq_genes) %>% 
  mutate(cumsum = cumsum(propability)) -> data_pie_chart


n <- 12
custom_palette <-c(brewer.pal(n = n + 1, name = "Paired")[-c(11, 12)], "white")

scales::show_col(ggsci::pal_npg("nrc", alpha = 1)(10))
# 
# scales::show_col(ggsci::pal_igv()(10))
# 
# custom_palette <- c(ggsci::pal_npg()(10), "white")

# ggsci::pal_npg("nrc")


svg("results/figures/gene-database-summary/filtered-1-12-gene-database-summary/piechart-filt-1-12-genes-occurence.svg", width = 8, height = 8)
ggplot(data_pie_chart, aes(x = "", y = n_genes, fill = freq_genes)) +
  geom_bar(stat = "identity", width = 1, color = "black", alpha = 0.6) +
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_manual(values = custom_palette) +
  labs(fill = "Gene Frequency") +
  ggtitle("Distribution of Gene Frequency")
dev.off()
