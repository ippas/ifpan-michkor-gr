# Przypisujemy wynik przetwarzania danych do obiektu df
df <- universal_gr_genes %>%
  mutate(label_plot = "universal") %>% 
  rbind(., universal_gr_genes %>% mutate(label_plot = simple_tissue)) %>% 
  ungroup() %>% 
  select(c("simple_tissue",  "regulation", "hgnc_symbol", 
           "label_regulation", "source", "log2ratio", "treatment", "dose", 
           "time", "treatment_type", "environment", "comparison", "n_genes", 
           "rank_score", "n_list", "n_genes_tissue", "sum_rs", "label_plot")) %>% 
  mutate(color = ifelse(regulation == "up", "firebrick", "blue4")) %>% 
  rbind(., top50_tissue_data %>%
          ungroup() %>%
          select(c("simple_tissue",  "regulation", "hgnc_symbol",
                   "label_regulation", "source", "log2ratio", "treatment", "dose",
                   "time", "treatment_type", "environment", "comparison", "n_genes",
                   "rank_score", "n_list", "n_genes_tissue", "sum_rs")) %>%
          mutate(label_plot = simple_tissue) %>%
          mutate(color = "gray80") %>% 
          filter(!(hgnc_symbol %in% universal_gr_genes$hgnc_symbol))
  ) %>%
  mutate(color = factor(color, levels = c("firebrick", "blue4", "gray80")),
         label_plot = factor(label_plot, levels = c("universal", "lung", "blood", 
                                                    "bone", "embryos", "kidney", 
                                                    "adrenal-gland", "adipose", 
                                                    "muscle", "brain", "liver", 
                                                    "cartilage",
                                                    "small-intestine", 
                                                    "placenta", "spleen")))


svg(filename = "results/figures/harmonized-tissue-gene-lists/summary-tissue-gene-lists/boxplot-universal-genes.svg", width = 14, height = 10)
# Tworzenie wykresu
ggplot(df, aes(x = label_plot, y = log2ratio, color = color)) +
  geom_jitter(alpha = 0.6, size = 1) +
  geom_boxplot(data = df %>% filter(label_plot == "universal"),
               aes(x = label_plot, y = log2ratio, group = regulation, fill = regulation),
               outlier.shape = NA, color = "black", size = 0.5,
               position = position_identity(), alpha = 0) +  # Boxplot tylko dla "universal"
  scale_fill_manual(values = c("up" = "firebrick", "down" = "blue4")) + # Kolory dla up i down
  scale_color_identity() +  # Zapewnia użycie kolorów z kolumny "color"
  theme_minimal() +
  theme(legend.position = "bottom") 
dev.off()


universal_gr_genes %>% 
  filter(simple_tissue == "spleen") %>% 
  select(c(regulation, simple_tissue, hgnc_symbol, sum_rs)) %>% 
  arrange(desc(sum_rs)) %>% 
  unique %>% as.data.frame()
