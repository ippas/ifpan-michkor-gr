HumanProteinAtlas_rnaTissueConsensus %>% 
  select(c(gene_symbol, Organ, organ_nTPM_sum, organ_nTPM_mean, organ_nTPM_median, organ_nTPM_max)) %>% 
  unique


HumanProteinAtlas_rnaTissueConsensus %>% 
  select(gene_symbol, Organ, organ_nTPM_sum, organ_nTPM_mean, organ_nTPM_median, organ_nTPM_max) %>% 
  distinct() %>% 
  pivot_longer(
    cols = starts_with("organ_nTPM"),
    names_to = "metric",
    values_to = "value"
  ) %>% 
  ggplot(aes(x = Organ, y = log2(value))) +
  geom_boxplot(outlier.alpha = 0.2) +
  facet_wrap(~ metric, scales = "free_y") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    x = "Organ",
    y = "nTPM",
    title = "Distribution of nTPM metrics across organs"
  )

gene_name <- "PER1"

HumanProteinAtlas_rnaTissueConsensus %>% 
  select(c(gene_symbol, Organ, organ_nTPM_sum, organ_nTPM_mean, organ_nTPM_median, organ_nTPM_max)) %>% 
  unique %>% 
  filter(gene_symbol == gene_name) %>% 
  select(Organ, organ_nTPM_sum, organ_nTPM_mean, organ_nTPM_median, organ_nTPM_max) %>% 
  pivot_longer(
    cols = starts_with("organ_nTPM"),
    names_to = "metric",
    values_to = "value"
  ) %>% 
  ggplot(aes(x = Organ, y = value)) +
  geom_col(
    fill = "darkred",
    color = "black",
    size = 1
  ) +
  facet_wrap(~ metric, scales = "free_y") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = paste("Expression of", gene_name, "across organs"),
    x = "Organ",
    y = "nTPM"
  )


cleanOverlapResults_grSystemicTissues_ICD10F %>% 
  # filter(grSignature == "bloodDown") %>%
  filter(grepl("blood", grSignature)) %>% 
  filter(p_value < 0.05) %>% 
  filter(observed_overlap >= 3) %>% .$overlap_genes %>% 
  strsplit(",") %>% 
  unlist %>% unique -> tmp


HumanProteinAtlas_rnaTissueConsensus %>% 
  filter(Organ %in% c("brain", "bone marrow & lymphoid tissues", "respiratory system")) %>% 
  select(gene_symbol, Organ, organ_nTPM_max) %>% 
  filter(gene_symbol %in% tmp) %>% 
  distinct() %>% 
  ggplot(aes(x = Organ, y = log2(organ_nTPM_max))) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.2) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    x = "Organ",
    y = "max nTPM",
    title = "Distribution of max nTPM across organs"
  )


df <- HumanProteinAtlas_rnaTissueConsensus %>% 
  filter(Organ %in% c("brain", "bone marrow & lymphoid tissues", "respiratory system")) %>% 
  select(gene_symbol, Organ, organ_nTPM_max) %>% 
  filter(gene_symbol %in% tmp) %>% 
  distinct() %>% 
  mutate(log_expr = log2(organ_nTPM_max + 1))

# one-way ANOVA
anova_res <- df %>%
  anova_test(log_expr ~ Organ)

anova_res


posthoc_res <- df %>%
  tukey_hsd(log_expr ~ Organ)

posthoc_res

anova_res <- df %>%
  anova_test(
    dv = log_expr,
    wid = gene_symbol,
    within = Organ
  )

anova_res


posthoc_res <- df %>%
  pairwise_t_test(
    log_expr ~ Organ,
    paired = TRUE,
    p.adjust.method = "BH"
  )

posthoc_res

