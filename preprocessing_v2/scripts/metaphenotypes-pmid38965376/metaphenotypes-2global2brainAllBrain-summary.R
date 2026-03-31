metaphenotypes_AllBrain2Brain2Global$overlap$original_data$df %>% 
  filter(!(Var1 %in% c("brain_down", "brain_up", "metasignature_up", "metasignature_down"))) %>% 
  select(-fdr) %>% 
  group_by(Var2) %>% 
  nest() %>% 
  .[1, 2] %>% .[[1]] %>% .[[1]] %>% 
  select(overlap_genes) %>% 
  filter(overlap_genes != "") %>% 
  .$overlap_genes %>% 
  strsplit(split = ",") %>% 
  unlist %>% unique() %>% length()
  
metaphenotypes_AllBrain2Brain2Global$overlap$original_data$df %>% 
  filter(!(Var1 %in% c("brain_down", "brain_up", "metasignature_up", "metasignature_down"))) %>% 
  select(-fdr) %>% 
  group_by(Var2) %>% 
  nest() %>% 
  mutate(
    n_unique_genes = map_int(data, ~ .x %>% 
                               filter(overlap_genes != "") %>% 
                               pull(overlap_genes) %>% 
                               strsplit(",") %>% 
                               unlist() %>% 
                               trimws() %>% 
                               unique() %>% 
                               length()
    )
  ) %>% 
  mutate(
    n_unique_genes_p0.05 = map_int(data, ~ .x %>% 
                               filter(p_value < 0.01) %>% 
                               filter(overlap_genes != "") %>% 
                               pull(overlap_genes) %>% 
                               strsplit(",") %>% 
                               unlist() %>% 
                               trimws() %>% 
                               unique() %>% 
                               length()
    )
  ) %>% 
  left_join(., data.frame(name_list = names(metaphenotypes_AllBrain2Brain2Global$overlap$gene_list_sizes),
                          n_phenotypes_genes = metaphenotypes_AllBrain2Brain2Global$overlap$gene_list_sizes %>% unname()),
            by = c("Var2" = "name_list")) -> metaphenotypes_AllBrain2Brain2Global$overlap$original_data$grouped_data 


metaphenotypes_AllBrain2Brain2Global$overlap$original_data$grouped_data %>% 
  mutate(no_signif_genes = n_unique_genes - n_unique_genes_p0.05) %>% 
  mutate(prop_n_signif_genes = n_unique_genes_p0.05/n_phenotypes_genes) %>% 
  mutate(prop_n_noSignif_genes = no_signif_genes/n_phenotypes_genes) %>% 
  mutate(prop_all_genes = n_unique_genes/n_phenotypes_genes)


df <- metaphenotypes_AllBrain2Brain2Global$overlap$original_data$grouped_data %>% 
  mutate(no_signif_genes = n_unique_genes - n_unique_genes_p0.05) %>% 
  mutate(prop_n_signif_genes = n_unique_genes_p0.05/n_phenotypes_genes) %>% 
  mutate(prop_n_noSignif_genes = no_signif_genes/n_phenotypes_genes) %>% 
  mutate(prop_all_genes = n_unique_genes/n_phenotypes_genes)

ggplot(df, aes(x = Var2, y = prop_n_signif_genes)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(
    x = "Var2",
    y = "Proportion of all genes",
    title = "Barplot of prop_all_genes per Var2"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

metaphenotypes_AllBrain2Brain2Global$overlap$original_data$grouped_data$Var2

metaphenotypes_AllBrain2Brain2Global$overlap$original_data$grouped_data %>% 
  filter(Var2 %in% c("alcohol_use_and_misuse", "anxiety_and_nervousness",
                     "cognition_and_processing_speed", "depressive_symptomatology",
                     "trauma")) %>% 
  select(-c(n_unique_genes, n_unique_genes_p0.05, n_phenotypes_genes)) %>% 
  unnest() %>% 
  filter(p_value < 0.05) %>% 
  as.data.frame() %>% 
  .$overlap_genes %>% 
  strsplit(",") %>% unlist %>% cat(sep = "\n")
