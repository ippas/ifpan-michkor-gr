GWASCatalogDisGeNETgenebass_GrTissue_overlapChi2$processed$original_data$df %>% 
  filter(grepl("_F", Var1)) %>% 
  filter(p_value < 0.05) %>% 
  filter(gene_overlap_count > 2) %>% 
  filter(grepl("_F4", Var1))


GWASCatalogDisGeNETgenebass_GrTissue_overlapChi2$processed$original_data$df %>%
  filter(grepl("_F", Var1)) %>%
  filter(p_value < 0.05) %>%
  filter(gene_overlap_count > 2) %>%
  filter(log2_odds_ratio > 0) %>% 
  mutate(
    F_category = factor(
      str_extract(Var1, "F[0-9]x"),
      levels = paste0("F", 0:9, "x")
    )
  ) %>%
  filter(Var1 != "GWASCatalog_F9x_F3x_F8x_F2x_attention deficit hyperactivity disorder,bipolar disorder,autism spectrum disorder,schizophrenia,major depressive disorder") %>% 
  filter(Var1 != "genebass_F9x_Date_F99_first_reported_(mental_disorder,_not_otherwise_specified)") %>% 
  na.omit() %>%   
  filter(grepl("F3", Var1)) %>%
  filter(grepl("lung", Var2, ignore.case = T)) %>% .$overlap_genes %>% unique %>% strsplit(",") %>% unlist %>% unique
  group_by(F_category) %>% 
  nest() %>% 
  mutate(
    n_association = map_int(data, nrow),
    n_genes = map_int(data, ~ .x$overlap_genes %>%
                        strsplit(",") %>%
                        unlist() %>%
                        unique() %>%
                        length())
  ) %>%
  select(-data) %>%
  ungroup() %>%
  
  ## <<< KLUCZOWE >>>
  complete(
    F_category = factor(paste0("F", 0:9, "x"), levels = paste0("F", 0:9, "x")),
    fill = list(n_association = 0, n_genes = 0)
  ) %>%
  arrange(F_category)



df_plot

brain_gene_lists <- papers_data_preprocessing %>% 
  filter(simple_tissue == "brain", n_genes >= 10) %>% 
  select(label, hgnc_symbol) %>% 
  group_by(label) %>% 
  summarise(genes = list(unique(hgnc_symbol)), .groups = "drop") %>% 
  deframe()



