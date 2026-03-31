# n categories
genebass_mentalHealth_skat %>% .$category %>% unique %>% length()

# n phenotypes
genebass_mentalHealth_skat %>% .$description %>% unique %>% length()

# n genes 
genebass_mentalHealth_skat$gene_symbol %>% unique %>% length()


genebass_mentalHealth_skat %>% 
  filter(gene_symbol %in% hgnc_symbols_vector_v110) %>% .$gene_symbol %>% 
  unique %>% length()


genebass_mentalHealth_skat %>% 
  filter(gene_symbol %in% hgnc_symbols_vector_v110) %>% 
  dim


genebass_mentalHealth_skat %>%
  group_by(category) %>%
  summarise(
    n_phenotypes = n_distinct(description),
    n_genes_total = n_distinct(gene_symbol),
    n_genes_in_biomart = n_distinct(gene_symbol[gene_symbol %in% hgnc_symbols_vector_v110]),
    n_rows_in_biomart = sum(gene_symbol %in% hgnc_symbols_vector_v110)
  ) %>%
  arrange(desc(n_phenotypes))


genebass_mentalHealth_skat %>% 
  filter(gene_symbol %in% hgnc_symbols_vector_v110) %>%
  group_by(category, description) %>%
  summarise(
    n_genes = n_distinct(gene_symbol),
    n_rsid = nrow(cur_data()),
    .groups = "drop"
  ) %>%
  nest(data = c(description, n_genes, n_rsid)) %>%
  mutate(
    n_phenotypes = map_int(data, nrow),
    n_genes_in_biomart = map_int(data, ~ sum(.x$n_genes)),
    n_rows_in_biomart = map_int(data, ~ sum(.x$n_rsid)),
    mean_genes_per_phenotype = map_dbl(data, ~ mean(.x$n_genes)),
    sd_genes_per_phenotype = map_dbl(data, ~ sd(.x$n_genes)),
    mean_rsid_per_phenotype = map_dbl(data, ~ mean(.x$n_rsid)),
    sd_rsid_per_phenotype = map_dbl(data, ~ sd(.x$n_rsid))
  ) %>%
  select(-data) %>%
  arrange(desc(n_phenotypes)) %>% 
  write_xlsx(
    path = "data/genebass/genebass_mental_health_summary_subcategories.xlsx"
  )
