AllBrainGeneDf

# n records
AllBrainGeneDf %>% nrow()

# n genes
AllBrainGeneDf %>% .$hgnc_symbol %>% unique() %>% length()

# n publications
AllBrainGeneDf$source %>% unique() %>% length()

# n geneLists
AllBrainGeneDf$label %>% unique() %>% length()

# summary freq of genes per list
AllBrainGeneDf$hgnc_symbol %>% 
  table() %>%
  as.data.frame() %>% 
  setNames(c("hgnc_symbol", "freq")) %>% 
  group_by(freq) %>% 
  nest() %>% 
  arrange(freq) %>% 
  mutate(
    n_genes = map_int(data, nrow),
    data = map(data, ~ .x %>% pull(hgnc_symbol) %>% as.character)   # 👈 zamiana tibble na wektor
  ) %>% 
  ungroup() %>% 
  mutate(
    fraction_genes = n_genes / sum(n_genes),
    cumulative_fraction = cumsum(fraction_genes)
  ) %>% 
  select(freq, n_genes, fraction_genes, cumulative_fraction, data)

# summary freq of genes per paper
AllBrainGeneDf %>%
  select(source, hgnc_symbol) %>% 
  unique() %>%
  .$hgnc_symbol %>% 
  table() %>%
  as.data.frame() %>% 
  setNames(c("hgnc_symbol", "freq")) %>% 
  group_by(freq) %>% 
  nest() %>% 
  arrange(freq) %>% 
  mutate(
    n_genes = map_int(data, nrow),
    data = map(data, ~ .x %>% pull(hgnc_symbol) %>% as.character)   # 👈 zamiana tibble na wektor
  ) %>% 
  ungroup() %>% 
  mutate(
    fraction_genes = n_genes / sum(n_genes),
    cumulative_fraction = cumsum(fraction_genes)
  ) %>% 
  select(freq, n_genes, fraction_genes, cumulative_fraction, data)

AllBrainGeneDf %>%
  select(regulation, hgnc_symbol) %>% 
  unique() %>% 
  group_by(hgnc_symbol) %>% 
  nest() %>% 
  mutate(n_regulation = map(data, ~ .x %>% nrow)) %>% 
  unnest(n_regulation) %>% 
  unnest(data) %>% 
  mutate(regulation = ifelse(
    n_regulation == 2,
    "both",
    ifelse(n_regulation == 1 &
             regulation == "up", "only_up", "only_down")
  )) %>% 
  select(-n_regulation) %>% 
  group_by(regulation) %>% 
  nest() %>% 
  mutate(n_genes = map(data, ~ .x %>% unique %>% nrow())) %>%
  mutate(genes = map(data, ~ .x$hgnc_symbol %>% unique %>% as.character)) %>% 
  unnest(n_genes) %>% 
  select(-data)
  
# ##############################################################################
# ---- summary upregulated genes ----
# ##############################################################################
# n records
AllBrainGeneDf %>% 
  filter(regulation == "up") %>% 
  nrow()

# n genes
AllBrainGeneDf %>% 
  filter(regulation == "up") %>% 
  .$hgnc_symbol %>% unique() %>% length()

# n publications
AllBrainGeneDf %>% 
  filter(regulation == "up") %>% 
  .$source %>% 
  unique() %>% length()

# n geneLists
AllBrainGeneDf %>% 
  filter(regulation == "up") %>% 
  .$label %>% 
  unique() %>% 
  length()


# summary freq of genes per list
AllBrainGeneDf %>% 
  filter(regulation == "up") %>% 
  .$hgnc_symbol %>% 
  table() %>%
  as.data.frame() %>% 
  setNames(c("hgnc_symbol", "freq")) %>% 
  group_by(freq) %>% 
  nest() %>% 
  arrange(freq) %>% 
  mutate(
    n_genes = map_int(data, nrow),
    data = map(data, ~ .x %>% pull(hgnc_symbol) %>% as.character)   # 👈 zamiana tibble na wektor
  ) %>% 
  ungroup() %>% 
  mutate(
    fraction_genes = n_genes / sum(n_genes),
    cumulative_fraction = cumsum(fraction_genes)
  ) %>% 
  select(freq, n_genes, fraction_genes, cumulative_fraction, data)

AllBrainGeneDf %>%
  filter(regulation == "up") %>% 
  select(source, hgnc_symbol) %>% 
  unique() %>%
  .$hgnc_symbol %>% 
  table() %>%
  as.data.frame() %>% 
  setNames(c("hgnc_symbol", "freq")) %>% 
  group_by(freq) %>% 
  nest() %>% 
  arrange(freq) %>% 
  mutate(
    n_genes = map_int(data, nrow),
    data = map(data, ~ .x %>% pull(hgnc_symbol) %>% as.character)   # 👈 zamiana tibble na wektor
  ) %>% 
  ungroup() %>% 
  mutate(
    fraction_genes = n_genes / sum(n_genes),
    cumulative_fraction = cumsum(fraction_genes)
  ) %>% 
  select(freq, n_genes, fraction_genes, cumulative_fraction, data)

# ##############################################################################
# ---- summary downregulated genes ----
# ##############################################################################
# n records
AllBrainGeneDf %>% 
  filter(regulation == "down") %>% 
  nrow()

# n genes
AllBrainGeneDf %>% 
  filter(regulation == "down") %>% 
  .$hgnc_symbol %>% unique() %>% length()

# n publications
AllBrainGeneDf %>% 
  filter(regulation == "down") %>% 
  .$source %>% 
  unique() %>% length()

# n geneLists
AllBrainGeneDf %>% 
  filter(regulation == "down") %>% 
  .$label %>% 
  unique() %>% 
  length()


# summary freq of genes per list
AllBrainGeneDf %>% 
  filter(regulation == "down") %>% 
  .$hgnc_symbol %>% 
  table() %>%
  as.data.frame() %>% 
  setNames(c("hgnc_symbol", "freq")) %>% 
  group_by(freq) %>% 
  nest() %>% 
  arrange(freq) %>% 
  mutate(
    n_genes = map_int(data, nrow),
    data = map(data, ~ .x %>% pull(hgnc_symbol) %>% as.character)   # 👈 zamiana tibble na wektor
  ) %>% 
  ungroup() %>% 
  mutate(
    fraction_genes = n_genes / sum(n_genes),
    cumulative_fraction = cumsum(fraction_genes)
  ) %>% 
  select(freq, n_genes, fraction_genes, cumulative_fraction, data)

AllBrainGeneDf %>%
  filter(regulation == "down") %>% 
  select(source, hgnc_symbol) %>% 
  unique() %>%
  .$hgnc_symbol %>% 
  table() %>%
  as.data.frame() %>% 
  setNames(c("hgnc_symbol", "freq")) %>% 
  group_by(freq) %>% 
  nest() %>% 
  arrange(freq) %>% 
  mutate(
    n_genes = map_int(data, nrow),
    data = map(data, ~ .x %>% pull(hgnc_symbol) %>% as.character)   # 👈 zamiana tibble na wektor
  ) %>% 
  ungroup() %>% 
  mutate(
    fraction_genes = n_genes / sum(n_genes),
    cumulative_fraction = cumsum(fraction_genes)
  ) %>% 
  select(freq, n_genes, fraction_genes, cumulative_fraction, data)



summarize_gene_dataset(df = AllBrainGeneDf)
