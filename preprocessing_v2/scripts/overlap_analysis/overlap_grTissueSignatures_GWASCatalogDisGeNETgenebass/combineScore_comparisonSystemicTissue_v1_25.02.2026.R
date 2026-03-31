grSystemic_GWASCatalogDisGeNETgenebass_list
grNeural_GWASCatalogDisGeNETgenebass_list
grBlood_GWASCatalogDisGeNETgenebass_list
grLung_GWASCatalogDisGeNETgenebass_list


bind_rows(grSystemic_GWASCatalogDisGeNETgenebass_list)
bind_rows(grNeural_GWASCatalogDisGeNETgenebass_list)
bind_rows(grBlood_GWASCatalogDisGeNETgenebass_list)
bind_rows(grLung_GWASCatalogDisGeNETgenebass_list)


grSystemicTissues_all <- bind_rows(
  bind_rows(grSystemic_GWASCatalogDisGeNETgenebass_list) %>% mutate(tissue = "Systemic"),
  bind_rows(grNeural_GWASCatalogDisGeNETgenebass_list)   %>% mutate(tissue = "Neural"),
  bind_rows(grBlood_GWASCatalogDisGeNETgenebass_list)    %>% mutate(tissue = "Blood"),
  bind_rows(grLung_GWASCatalogDisGeNETgenebass_list)     %>% mutate(tissue = "Lung")
)


# ---- systemic ----
grSystemicTissues_all %>% 
  filter(p_value < 0.05) %>% 
  filter(observed_overlap >= 3) %>% 
  mutate(label = paste0(phenotype, "_", source)) %>% 
  filter(grSignature %in% c("systemicUp", "systemicDown")) %>% 
  .$label %>% unique -> all_signif_phenotype


grSystemicTissues_all %>% 
  mutate(label = paste0(phenotype, "_", source)) %>% 
  filter(label %in% all_signif_phenotype,
         grSignature %in% c("systemicUp", "systemicDown")) %>% 
  ggplot(aes(x = grSignature, y = combine_score)) +
  geom_boxplot(alpha = 0.4, width = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = grSignature),
              width = 0.12,
              size = 2,
              alpha = 0.7) +
  scale_color_manual(values = c("systemicUp" = "firebrick3",
                                "systemicDown" = "navy")) +
  theme_classic() +
  theme(legend.position = "none")

grSystemicTissues_all %>% 
  mutate(label = paste0(phenotype, "_", source)) %>% 
  filter(label %in% all_signif_phenotype) %>% 
  filter(grSignature %in% c("systemicUp", "systemicDown")) %>% 
  select(label, grSignature, combine_score) %>% 
  pivot_wider(names_from = grSignature,
              values_from = combine_score) %>% 
  drop_na() %>% 
  with(wilcox.test(systemicUp, systemicDown, paired = TRUE))

# ---- neural ----
grSystemicTissues_all %>% 
  filter(p_value < 0.05) %>% 
  filter(observed_overlap >= 3) %>% 
  mutate(label = paste0(phenotype, "_", source)) %>% 
  filter(grSignature %in% c("neuralUp", "neuralDown")) %>% 
  .$label %>% unique -> all_signif_phenotype


grSystemicTissues_all %>% 
  mutate(label = paste0(phenotype, "_", source)) %>% 
  filter(label %in% all_signif_phenotype,
         grSignature %in% c("neuralUp", "neuralDown")) %>% 
  ggplot(aes(x = grSignature, y = combine_score)) +
  geom_boxplot(alpha = 0.4, width = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = grSignature),
              width = 0.12,
              size = 2,
              alpha = 0.7) +
  scale_color_manual(values = c("neuralUp" = "firebrick3",
                                "neuralDown" = "navy")) +
  theme_classic() +
  theme(legend.position = "none")

grSystemicTissues_all %>% 
  mutate(label = paste0(phenotype, "_", source)) %>% 
  filter(label %in% all_signif_phenotype) %>% 
  filter(grSignature %in% c("neuralUp", "neuralDown")) %>% 
  select(label, grSignature, combine_score) %>% 
  pivot_wider(names_from = grSignature,
              values_from = combine_score) %>% 
  drop_na() %>% 
  with(wilcox.test(neuralUp, neuralDown, paired = TRUE))


# ---- blood ----
grSystemicTissues_all %>% 
  filter(p_value < 0.05) %>% 
  filter(observed_overlap >= 3) %>% 
  mutate(label = paste0(phenotype, "_", source)) %>% 
  filter(grSignature %in% c("bloodUp", "bloodDown")) %>% 
  .$label %>% unique -> all_signif_phenotype


grSystemicTissues_all %>% 
  mutate(label = paste0(phenotype, "_", source)) %>% 
  filter(label %in% all_signif_phenotype,
         grSignature %in% c("bloodUp", "bloodDown")) %>% 
  ggplot(aes(x = grSignature, y = combine_score)) +
  geom_boxplot(alpha = 0.4, width = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = grSignature),
              width = 0.12,
              size = 2,
              alpha = 0.7) +
  scale_color_manual(values = c("bloodUp" = "firebrick3",
                                "bloodDown" = "navy")) +
  theme_classic() +
  theme(legend.position = "none")

grSystemicTissues_all %>% 
  mutate(label = paste0(phenotype, "_", source)) %>% 
  filter(label %in% all_signif_phenotype) %>% 
  filter(grSignature %in% c("bloodUp", "bloodDown")) %>% 
  select(label, grSignature, combine_score) %>% 
  pivot_wider(names_from = grSignature,
              values_from = combine_score) %>% 
  drop_na() %>% 
  with(wilcox.test(bloodUp, bloodDown, paired = TRUE))

# ---- lung ----
grSystemicTissues_all %>% 
  filter(p_value < 0.05) %>% 
  filter(observed_overlap >= 3) %>% 
  mutate(label = paste0(phenotype, "_", source)) %>% 
  filter(grSignature %in% c("lungUp", "lungDown")) %>% 
  .$label %>% unique -> all_signif_phenotype


grSystemicTissues_all %>% 
  mutate(label = paste0(phenotype, "_", source)) %>% 
  filter(label %in% all_signif_phenotype,
         grSignature %in% c("lungUp", "lungDown")) %>% 
  ggplot(aes(x = grSignature, y = combine_score)) +
  geom_boxplot(alpha = 0.4, width = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = grSignature),
              width = 0.12,
              size = 2,
              alpha = 0.7) +
  scale_color_manual(values = c("lungUp" = "firebrick3",
                                "lungDown" = "navy")) +
  theme_classic() +
  theme(legend.position = "none")

grSystemicTissues_all %>% 
  mutate(label = paste0(phenotype, "_", source)) %>% 
  filter(label %in% all_signif_phenotype) %>% 
  filter(grSignature %in% c("lungUp", "lungDown")) %>% 
  select(label, grSignature, combine_score) %>% 
  pivot_wider(names_from = grSignature,
              values_from = combine_score) %>% 
  drop_na() %>% 
  with(wilcox.test(lungUp, lungDown, paired = TRUE))
  
