plot_gene_expression_boxplot <- function(data, hgnc_symbols) {
  data %>% 
    filter(regulation %in% c("down", "up")) %>% 
    filter(!log2ratio %in% c("", "NA", "Inf", "-Inf")) %>%
    filter(hgnc_symbol %in% hgnc_symbols) %>% 
    filter(!is.na(hgnc_symbol)) %>% 
    mutate(
      log2ratio = as.numeric(log2ratio),
      hgnc_symbol = fct_reorder(hgnc_symbol, hgnc_occurence, .desc = TRUE)
    ) %>% 
    ggplot(aes(x = hgnc_label, y = log2ratio)) +
    geom_boxplot(outlier.shape = NA, size = 1) +
    geom_jitter(width = 0.2, alpha = 1, color = "black") +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "firebrick", size = 1) +
    geom_hline(yintercept = -0.5, linetype = "dashed", color = "blue4", size = 1) +
    
    # Dodanie etykiet dla linii referencyjnych na początku osi X z log2 jako indeks dolny
    annotate("text", x = -Inf, y = 0.5, 
             label = expression(log[2] ~ "(FC) = 0.5"), 
             color = "firebrick", size = 5, hjust = 0, vjust = -0.5) +
    annotate("text", x = -Inf, y = -0.5, 
             label = expression(log[2] ~ "(FC) = -0.5"), 
             color = "blue4", size = 5, hjust = 0, vjust = 1.5) +
    
    theme_minimal() +
    labs(
      x = NULL,  # Usunięcie tytułu osi X
      y = "log2ratio"
    ) +
    theme(
      axis.text.x = element_text(size = 18, lineheight = 1, color = "black"),
      axis.text.y = element_text(size = 16, color = "black"),
      axis.title.y = element_text(size = 18, margin = margin(r = 10))  # Zwiększona odległość tytułu Y
    ) +
    coord_cartesian(clip = "off")
}

plot_gene_expression_boxplot(
  data = filt_1_12_gr_database,
  hgnc_symbols = c("FKBP5", "TSC22D3", "KLF9", "PER1", "PDK4")
)



filt_1_12_gr_database %>% 
  filter(regulation %in% c("down", "up")) %>% 
  filter(!log2ratio %in% c("", "NA", "Inf", "-Inf")) %>%
  filter(!is.na(log2ratio)) %>% 
  # filter(hgnc_symbol %in% c("FKBP5", "TSC22D3", "KLF9", "PER1", "PDK4")) %>% 
  group_by(hgnc_symbol, hgnc_occurence) %>% 
  nest() %>% 
  mutate(log2ratio_hgnc_occurence = map(data, ~nrow(.x))) %>% 
  unnest(log2ratio_hgnc_occurence) %>% 
  # filter(log2ratio_hgnc_occurence > 20) %>% 
  mutate(variance_expression = map(data, ~var(.x$log2ratio))) %>% 
  mutate(mean_expression = map_dbl(data, ~ mean(as.numeric(as.character(.x$log2ratio)), na.rm = TRUE))) %>% 
  unnest(variance_expression) %>% 
  unnest(mean_expression) %>% 
  mutate(cv = variance_expression/mean_expression) %>% 
  arrange(desc(cv)) %>% 
  mutate(n_down = map(data, ~nrow(filter(.x, regulation == "down")))) %>% 
  mutate(n_up = map(data, ~nrow(filter(.x, regulation == "up")))) %>% 
  unnest(c(n_down, n_up)) %>% 
  arrange(desc(n_down)) %>% 
  mutate(up_perc = n_up/log2ratio_hgnc_occurence) %>% 
  mutate(down_perc = n_down/log2ratio_hgnc_occurence) -> tmp

tmp %>%
  mutate(q1 = map(data, ~ quantile(as.numeric(.x$log2ratio), probs = 0.25, na.rm = TRUE))) %>% 
  mutate(q3 = map(data, ~ quantile(as.numeric(.x$log2ratio), probs = 0.75, na.rm = TRUE))) %>% 
  unnest(q1, q3) -> tmp


tmp %>% 
  mutate(
    new_regulation_perc = case_when(
      up_perc > 0.90 ~ "up",
      up_perc < 0.1 ~ "down",
      up_perc <= 0.90 & up_perc > 0.7 ~ "up_weak",
      up_perc >= 0.1 & up_perc < 0.3 ~ "down_weak",
      TRUE ~ "both"
    )
  ) %>% 
  mutate(
    new_regulation_quartile = case_when(
      q1 > 0.5 ~ "up",
      q3 < -0.5 ~ "down",
      q1 <= 0.5 & q1 > -0.5 ~ "up_weak",
      q3 >= -0.05 & q3 < 0.5 ~ "down_weak",
      TRUE ~ "both"
    )
  ) -> hgnc_new_regulation
  
hgnc_new_regulation %>%
  select(-data) %>% 
  as.data.frame() %>% 
  group_by(new_regulation_perc) %>% 
  slice_max(order_by = log2ratio_hgnc_occurence, n = 5) %>%
  filter(hgnc_symbol != "GEM") %>%
  ungroup %>% 
  select(hgnc_symbol, new_regulation_perc) %>% 
  split(.$new_regulation_perc) -> top5_hgnc_regulation_perc

hgnc_new_regulation %>%
  select(-data) %>% 
  as.data.frame() %>% 
  group_by(new_regulation_quartile) %>% 
  slice_max(order_by = log2ratio_hgnc_occurence, n = 5) %>%
  filter(hgnc_symbol != "NEDD9") %>%
  ungroup %>% 
  select(hgnc_symbol, new_regulation_quartile) %>% 
  split(.$new_regulation_quartile) -> top5_hgnc_regulation_quartile
  
# Określenie preferowanej kolejności
desired_order <- c("up", "up_weak", "both", "down_weak", "down")

# Przestawienie elementów listy według określonej kolejności
top5_hgnc_regulation_perc <- top5_hgnc_regulation_perc[desired_order]
top5_hgnc_regulation_quartile <- top5_hgnc_regulation_quartile[desired_order]

lapply(top5_hgnc_regulation_perc, function(x){
  plot_gene_expression_boxplot(
    data = filt_1_12_gr_database,
    hgnc_symbols = x$hgnc_symbol
  )
})
  

library(patchwork)

# Generowanie listy wykresów
plots <- lapply(top5_hgnc_regulation_perc, function(x) {
  plot_gene_expression_boxplot(
    data = filt_1_12_gr_database,
    hgnc_symbols = x$hgnc_symbol
  )
})

# Łączenie wykresów w jedną kolumnę
combined_plot <- Reduce(`+`, plots) + plot_layout(ncol = 1)

# Wyświetlenie wykresu
print(combined_plot)


# Generowanie listy wykresów
plots <- lapply(top5_hgnc_regulation_quartile, function(x) {
  plot_gene_expression_boxplot(
    data = filt_1_12_gr_database,
    hgnc_symbols = x$hgnc_symbol
  )
})

# Łączenie wykresów w jedną kolumnę
combined_plot <- Reduce(`+`, plots) + plot_layout(ncol = 1)

# Wyświetlenie wykresu
print(combined_plot)



# summary regulation perc
hgnc_new_regulation %>% 
  select(-c(log2ratio_hgnc_occurence, variance_expression, mean_expression, cv)) %>% 
  unnest(data) %>% 
  group_by(new_regulation_quartile) %>% 
  nest() %>% 
  mutate(n_hgnc = map(data, ~length(unique(.x$hgnc_symbol)))) %>% 
  mutate(n_papers = map(data, ~length(unique(.x$source)))) %>% 
  mutate(n_experiments = map(data, ~length(unique(.x$label)))) %>%
  mutate(n_tissues = map(data, ~length(unique(.x$simple_tissue)))) %>%
  mutate(mean = map(data, ~mean(as.numeric(.x$log2ratio)))) %>%
  mutate(median = map(data, ~median(as.numeric(.x$log2ratio)))) %>%
  mutate(variance = map(data, ~var(as.numeric(.x$log2ratio)))) %>%
  unnest(n_hgnc, n_papers, n_experiments, n_tissues, mean, median, variance) %>% 
  mutate(cv = variance/mean) %>%
  select(-data) %>% 
  head %>% 
  mutate(
    new_regulation_quartile = case_when(
      new_regulation_quartile == "up" ~ "strongly upregulated",
      new_regulation_quartile == "down" ~ "strongly downregulated",
      new_regulation_quartile == "up_weak" ~ "moderately upregulated",
      new_regulation_quartile == "down_weak" ~ "moderately downregulated",
      new_regulation_quartile == "both" ~ "bidirectionally regulated"
    )
  ) %>% 
  column_to_rownames(var = "new_regulation_quartile") %>% 
  t() %>% 
  as.data.frame() %>% 
  select("strongly upregulated", "moderately upregulated", "bidirectionally regulated", "moderately downregulated", "strongly downregulated")
  
hgnc_new_regulation %>% 
  select(-c(log2ratio_hgnc_occurence, variance_expression, mean_expression, cv)) %>% 
  unnest(data) %>% 
  group_by(new_regulation_perc) %>% 
  nest() %>% 
  mutate(n_hgnc = map(data, ~length(unique(.x$hgnc_symbol)))) %>% 
  mutate(n_papers = map(data, ~length(unique(.x$source)))) %>% 
  mutate(n_experiments = map(data, ~length(unique(.x$label)))) %>%
  mutate(n_tissues = map(data, ~length(unique(.x$simple_tissue)))) %>%
  mutate(mean = map(data, ~mean(as.numeric(.x$log2ratio)))) %>%
  mutate(median = map(data, ~median(as.numeric(.x$log2ratio)))) %>%
  mutate(variance = map(data, ~var(as.numeric(.x$log2ratio)))) %>%
  unnest(n_hgnc, n_papers, n_experiments, n_tissues, mean, median, variance) %>% 
  mutate(cv = variance/mean) %>%
  select(-data) %>% 
  head %>% 
  mutate(
    new_regulation_perc = case_when(
      new_regulation_perc == "up" ~ "strongly upregulated",
      new_regulation_perc == "down" ~ "strongly downregulated",
      new_regulation_perc == "up_weak" ~ "moderately upregulated",
      new_regulation_perc == "down_weak" ~ "moderately downregulated",
      new_regulation_perc == "both" ~ "bidirectionally regulated"
    )
  ) %>% 
  column_to_rownames(var = "new_regulation_perc") %>% 
  t() %>% 
  as.data.frame() %>% 
  select("strongly upregulated", "moderately upregulated", "bidirectionally regulated", "moderately downregulated", "strongly downregulated")


hgnc_new_regulation %>% 
  select(-c(log2ratio_hgnc_occurence, variance_expression, mean_expression, cv)) %>% 
  unnest(data) %>% 
  group_by(new_regulation_perc, new_regulation_quartile) %>% 
  nest() %>% 
  mutate(n_both = map(data, ~length(unique(.x$hgnc_symbol)))) %>% 
  unnest(n_both) %>% 
  filter(new_regulation_perc == new_regulation_quartile)

# n records down and up
filt_1_12_gr_database %>% 
  group_by(hgnc_symbol) %>% 
  nest() %>% 
  mutate(n_regulation = map(data, ~length(unique(.x$regulation)))) %>% 
  unnest(n_regulation) %>% 
  mutate(n_up = map(data, ~nrow(filter(.x, regulation == "up")))) %>% 
  mutate(n_down = map(data, ~nrow(filter(.x, regulation == "down")))) %>% 
  unnest(n_up, n_down) -> assess_regulation_records
  
assess_regulation_records %>% 
  filter(n_regulation == 1) %>% 
  .$n_up %>% sum
  
assess_regulation_records %>% 
  filter(n_regulation == 1) %>% 
  .$n_down %>% sum

assess_regulation_records %>% 
  filter(n_regulation == 2) %>% 
  .$n_up %>% sum

filt_1_12_gr_database %>% 
  filter(hgnc_symbol == "NET1") %>% 
  arrange(regulation)

plot_log2ratio_gr_database(
  data = filt_1_12_gr_database, 
  hgnc_symbol = c("PTGS2", "PER2", "CDKN1C", "DUSP5", "NET1"), 
  point_size = 2,
  order_by_gene = "PTGS2" 
)
