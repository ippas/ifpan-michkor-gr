overlap_raw_signatures$
overlap_minus_clusters$plot_or
overlap_minus_global6$significant_uniq_data


t.test(
  overlap_raw_signatures$processed$original_data$df %>% 
    remove_duplicate_pairs() %>% 
    .$odds_ratio,
  overlap_minus_global6$processed$original_data$df %>% 
    remove_duplicate_pairs() %>% 
    .$odds_ratio,
  var.equal = T
)



overlap_minus_global6$processed$original_data$df %>% 
  remove_duplicate_pairs() %>% 
  .$odds_ratio


overlap_minus_global6$plot_or
overlap_minus_clusters

bind_rows(
  overlap_raw_signatures$df %>% 
    remove_duplicate_pairs() %>% 
    mutate(group = "Raw"),
  overlap_minus_clusters$df %>% 
    remove_duplicate_pairs() %>% 
    mutate(group = "Minus genes from clusters K/P/O/D"),
  overlap_minus_global6$df %>% 
    remove_duplicate_pairs() %>% 
    mutate(group = "Minus global (6 tissues)"),
  overlap_minus_global5$df %>% 
    remove_duplicate_pairs() %>% 
    mutate(group = "Minus global (5 tissues)")
) %>%
  mutate(group = factor(group, levels = c(
    "Raw",
    "Minus genes from clusters K/P/O/D",
    "Minus global (6 tissues)",
    "Minus global (5 tissues)"
  ))) %>%
  ggplot(aes(x = group, y = OR)) +
  geom_boxplot(fill = "gray90", color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.5) +
  labs(
    x = NULL,
    y = "Odds Ratio (OR)",
    title = "Comparison of OR distributions across GR-signature datasets"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 20, hjust = 1)
  )

overlap_raw_signatures$df %>% 
  remove_duplicate_pairs()



# overlap_raw_signatures$processed$original_data$list$
