# Opis:
# Skrypt generuje eksploracyjne wykresy ilustrujące powiązania pomiędzy
# sygnaturami genów zależnych od receptora glukokortykoidowego (GR) a
# ukrytymi czynnikami fenotypowymi (latent factors).
#
# Zakres:
# - Barploty liczby genów i czynników
# - Boxploty p-wartości i efektów beta
# - Scatterploty (beta vs -log10(p))
# - Heatmapy asocjacji gen × faktor
# - Odległość SNP od TSS (log10)
#

grSignatureLite_association_factorsP1e4Tss100kb <- factors_rsidGenes_P1e4tss100kb %>% 
  inner_join(lite_grSignatures, by = "gene_symbol") %>%
  dplyr::select(signature_name, everything()) %>% 
  group_by(signature_name, factor_id, gene_symbol) %>%
  slice_min(pvalue, with_ties = FALSE) %>% 
  ungroup %>% 
  mutate(
    in_GWASCatalog = case_when(
      signature_name == "metasignature_up"   & gene_symbol %in% metasignatureUp_GWASCatalog   ~ 1,
      signature_name == "metasignature_down" & gene_symbol %in% metasignatureDown_GWASCatalog ~ 1,
      signature_name == "brain_up"           & gene_symbol %in% brainUp_GWASCatalog           ~ 1,
      signature_name == "brain_down"         & gene_symbol %in% brainDown_GWASCatalog         ~ 1,
      TRUE ~ 0
    )
  ) %>% 
  select(signature_name, gene_symbol, in_GWASCatalog, everything()) %>%
  arrange(pvalue) -> grSignatureLite_association_factorsP1e4Tss100kb_df

grSignatureLite_association_factorsP1e4Tss100kb_df %>% head

factors_rsidGenes_P1e4tss100kb %>%
  distinct(factor_id, gene_symbol) %>%
  count(factor_id, name = "n_genes") %>%
  ggplot(aes(x = factor_id, y = n_genes)) +
  geom_col(fill = "steelblue") +
  labs(
    title = "Number of GR-Dependent Genes per Factor",
    x = "Factor ID",
    y = "Number of Genes"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )


grSignatureLite_association_factorsP1e4Tss100kb_df %>%
  group_by(signature_name) %>%
  summarise(n_factors = n_distinct(factor_id)) %>%
  ungroup() %>% 
  ggplot(aes(x = signature_name, y = n_factors)) +
  geom_col(fill = "steelblue") +
  labs(
    title = "Number of Unique Factors per GR Signature",
    x = "GR Signature",
    y = "Number of Factors"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# Liczba unikalnych genów per sygnatura
grSignatureLite_association_factorsP1e4Tss100kb_df %>% 
  group_by(signature_name) %>%
  summarise(n_genes = n_distinct(gene_symbol)) %>%
  ungroup() %>% 
  ggplot(aes(x = signature_name, y = n_genes)) +
  geom_col(fill = "gray") +
  labs(
    title = "Number of Unique GR-Dependent Genes per Signature",
    x = "GR Signature",
    y = "Number of Genes"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Stałe poziomy faktorów
factor_levels <- paste0("factor_", 1:35)

# Wykres bez zmiennych pośrednich

grSignatureLite_association_factorsP1e4Tss100kb_df %>%
  distinct(signature_name, factor_id, gene_symbol) %>% 
  count(signature_name, factor_id, name = "n_genes") %>% 
  complete(signature_name, factor_id = factor_levels, fill = list(n_genes = 0)) %>%
  # mutate(factor_id = factor(factor_id, levels = factor_levels)) %>%.$factor_id %>% unique
  # filter(is.na(factor_id))
  ggplot(aes(x = factor_id, y = n_genes)) +
  geom_col(fill = "gray") +
  facet_wrap(~ signature_name, ncol = 1) +
  labs(
    title = "Number of GR-Dependent Genes per Factor",
    x = "Factor ID",
    y = "Number of Genes"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

grSignatureLite_association_factorsP1e4Tss100kb_df %>%
  filter(!is.na(pvalue), pvalue > 0) %>%
  mutate(log10_pvalue = -log10(pvalue)) %>%
  ggplot(aes(x = signature_name, y = log10_pvalue)) +
  geom_boxplot(outlier.shape = NA, fill = "skyblue") +
  geom_jitter(width = 0.2, alpha = 0.2, size = 0.8) +
  geom_text(
    data = function(d) d %>%
      filter(!is.na(pvalue), pvalue > 0, -log10(pvalue) > 25) %>%
      mutate(log10_pvalue = -log10(pvalue),
             label = paste0(gene_symbol, ", ", factor_id)),
    aes(label = label),
    size = 2.5,
    vjust = -0.3,
    check_overlap = TRUE
  ) +
  labs(
    title = "Distribution of -log10(p-value) per GR Signature",
    x = "GR Signature",
    y = "-log10(p-value)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


grSignatureLite_association_factorsP1e4Tss100kb_df %>%
  filter(!is.na(beta), !is.na(pvalue), pvalue > 0) %>%
  mutate(log10_pvalue = -log10(pvalue)) %>%
  ggplot(aes(x = beta, y = log10_pvalue)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_text(
    data = function(d) d %>%
      filter(!is.na(beta), !is.na(pvalue), pvalue > 0) %>%
      mutate(log10_pvalue = -log10(pvalue)) %>%
      filter(log10_pvalue > 20 | beta > 0.15 | beta < -0.15) %>%
      mutate(label = paste0(gene_symbol, ", ", factor_id)),
    aes(label = label),
    size = 2.5,
    vjust = -0.5,
    check_overlap = TRUE
  ) +
  facet_wrap(~ signature_name, ncol = 2) +
  labs(
    title = "Effect Size vs. -log10(p-value) per GR Signature",
    x = "Effect Size (Beta)",
    y = "-log10(p-value)"
  ) +
  theme_minimal()


grSignatureLite_association_factorsP1e4Tss100kb_df %>%
  filter(signature_name == "brain_down") %>%
  select(gene_symbol, factor_id, beta) %>%
  distinct() %>%
  pivot_wider(names_from = factor_id, values_from = beta) %>%
  column_to_rownames("gene_symbol") %>%
  as.matrix() %>%
  scale(center = TRUE, scale = FALSE) %>%  # opcjonalnie: centrowanie po genach
  {
    heatmap_data <- .
    ggplot(as.data.frame(as.table(heatmap_data)), aes(Var2, Var1, fill = Freq)) +
      geom_tile() +
      scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
      labs(
        title = "Heatmap of Beta Values (brain_down)",
        x = "Factor ID",
        y = "Gene Symbol",
        fill = "Beta"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size = 6)
      )
  }



grSignatureLite_association_factorsP1e4Tss100kb_df %>%
  select(signature_name, gene_symbol, factor_id, beta) %>%
  distinct() %>%
  group_by(signature_name) %>%
  complete(gene_symbol, factor_id) %>%  # teraz tylko w obrębie jednej sygnatury
  ungroup() %>%
  ggplot(aes(x = factor_id, y = gene_symbol, fill = beta)) +
  geom_tile(color = "grey90", linewidth = 0.1) +
  facet_wrap(~ signature_name, ncol = 2, scales = "free_y") +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    na.value = "grey90",  # kolor dla braków
    name = "Beta"
  ) +
  labs(
    title = "Heatmap of GR-Dependent Gene Associations per Signature",
    x = "Factor ID",
    y = "Gene Symbol"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 5),
    strip.text = element_text(face = "bold")
  )


grSignatureLite_association_factorsP1e4Tss100kb_df %>%
  filter(!is.na(distance_to_tss)) %>%
  mutate(
    abs_distance = abs(distance_to_tss),
    log10_distance = log10(abs_distance + 1)
  ) %>%
  ggplot(aes(x = log10_distance)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white") +
  facet_wrap(~ signature_name, ncol = 2) +
  labs(
    title = "Distribution of SNP Distance to TSS (log10 scale)",
    x = "log10(Distance to TSS + 1) [bp]",
    y = "Number of SNP-Gene Pairs"
  ) +
  theme_minimal()



grSignatureLite_association_factorsP1e4Tss100kb_df %>%
  filter(!is.na(distance_to_tss)) %>%
  mutate(
    abs_distance = abs(distance_to_tss),
    log10_distance = log10(abs_distance + 1)
  ) %>%
  ggplot(aes(x = "", y = log10_distance)) +
  geom_boxplot(fill = "darkorange", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.2, size = 0.8) +
  geom_text(
    data = function(d) d %>%
      filter(!is.na(distance_to_tss), abs(distance_to_tss) < 1000) %>%
      mutate(
        abs_distance = abs(distance_to_tss),
        log10_distance = log10(abs_distance + 1),
        label = paste0(gene_symbol, ", ", factor_id)
      ),
    aes(label = label),
    size = 2.5,
    vjust = -0.3,
    check_overlap = TRUE
  ) +
  facet_wrap(~ signature_name, ncol = 2) +
  labs(
    title = "Distance to TSS per GR Signature (with genes near TSS < 1kb)",
    x = NULL,
    y = "log10(Distance to TSS + 1) [bp]"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(face = "bold")
  )


grSignatureLite_association_factorsP1e4Tss100kb_df %>%
  filter(!is.na(distance_to_tss)) %>%
  mutate(
    abs_distance = abs(distance_to_tss),
    log10_distance = log10(abs_distance + 1)
  ) %>%
  ggplot(aes(x = "", y = log10_distance)) +
  geom_boxplot(fill = "darkorange", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.2, size = 0.8) +
  geom_text(
    data = function(d) d %>%
      filter(!is.na(distance_to_tss), abs(distance_to_tss) < 10000) %>%
      mutate(
        abs_distance = abs(distance_to_tss),
        log10_distance = log10(abs_distance + 1),
        label = paste0(gene_symbol, ", ", factor_id)
      ),
    aes(label = label),
    size = 2.5,
    vjust = -0.3,
    check_overlap = TRUE
  ) +
  facet_wrap(~ signature_name, ncol = 2) +
  labs(
    title = "Distance to TSS per GR Signature (genes <10kb from TSS)",
    x = NULL,
    y = "log10(Distance to TSS + 1) [bp]"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(face = "bold")
  )


