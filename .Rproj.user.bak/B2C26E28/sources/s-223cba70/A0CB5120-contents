################################################################################
# boxplot for tissue gene lists (top 50 per tissue and regulation) 
plot_data_sorted <- plot_data %>% 
  unnest(data) %>% 
  group_by(regulation, simple_tissue, hgnc_symbol) %>% 
  nest() %>% 
  mutate(sum_rs = map(data, ~ .x$rank_score %>% sum)) %>% 
  unnest(sum_rs) %>% 
  group_by(regulation, simple_tissue) %>% 
  nest() %>% 
  mutate(top_50_rs = map(data, ~ .x %>% arrange(desc(sum_rs)) %>% head(50))) %>% 
  select(-data) %>% 
  unnest(top_50_rs) %>% 
  unnest(data) %>%
  ungroup() %>% 
  group_by(simple_tissue) %>% 
  nest() %>% 
  mutate(n_records = map(data, ~nrow(.x))) %>% 
  unnest(n_records) %>% 
  arrange(desc(n_records)) 

# **Ustawienie poziomów faktora zgodnie z kolejnością posortowaną**
plot_data_sorted <- plot_data_sorted %>%
  mutate(simple_tissue = factor(simple_tissue, levels = plot_data_sorted$simple_tissue)) %>% unnest(data)

plot_data_sorted %>% ungroup %>% class

plot_data_sorted <- plot_data_sorted %>% 
  ungroup %>% 
  mutate(label = str_replace_all(label_regulation, "_down", "")) %>% 
  mutate(label = str_replace_all(label, "_up", "")) %>% 
  group_by(simple_tissue) %>% 
  nest %>% 
  mutate(n_expreriments = map(data, ~ .x$label %>% unique() %>% length)) %>% 
  mutate(n_lists = map(data, ~ .x$label_regulation %>% unique() %>% length)) %>% 
  mutate(n_genes_up = map(data, ~ .x %>% filter(regulation == "up") %>% .$hgnc_symbol %>% unique %>% length)) %>% 
  mutate(n_genes_down = map(data, ~ .x %>% filter(regulation == "down") %>% .$hgnc_symbol %>% unique %>% length)) %>% 
  mutate(n_records_up = map(data, ~ .x %>% filter(regulation == "up") %>% nrow())) %>% 
  mutate(n_records_down = map(data, ~ .x %>% filter(regulation == "down") %>% nrow())) %>% 
  unnest() %>% 
  ungroup() %>% 
  mutate(label_figure = paste(simple_tissue, n_expreriments, n_lists, n_genes_up, n_genes_down, n_records_up, n_records_down, sep = "\n")) 
# mutate(label_figure = factor(label_figure, levels = unique(label_figure[match(simple_tissue, levels(simple_tissue))])))

plot_data_sorted %>%
  group_by(simple_tissue, label_figure) %>% 
  nest() %>% 
  mutate(label_figure = factor(label_figure, levels = label_figure[order(match(simple_tissue, levels(simple_tissue)))])) %>% 
  unnest(data) -> plot_data_sorted

plot_data_sorted -> top50_tissue_data

ggplot(plot_data_sorted, aes(x = label_figure, y = log2ratio, color = regulation)) +
  geom_jitter(alpha = 0.6, size = 0.5) +
  geom_boxplot(aes(group = interaction(simple_tissue, regulation)), 
               outlier.shape = NA, color = "black", fill = NA, size = 0.5,
               position = "identity") +  # Boxploty bez przesunięcia
  theme_minimal() +
  theme(
    legend.position = "bottom",  # Przeniesienie legendy pod spód
    axis.text.x = element_text(hjust = 0, vjust = 0.5), # Rotacja tekstu
    axis.title.x = element_blank(),  # Usunięcie tytułu osi X, jeśli niepotrzebny
    axis.ticks.x.top = element_line(),  # Dodanie znaczników na górnej osi X
    axis.text.x.top = element_text(size = rel(1))  # Stylizacja górnej osi
  ) +
  scale_x_discrete(position = "top") +  # Przeniesienie osi X na górę
  scale_color_manual(values = c("down" = "blue4", "up" = "firebrick"))
