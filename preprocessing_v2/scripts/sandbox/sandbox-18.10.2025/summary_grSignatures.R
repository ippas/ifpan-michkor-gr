flat_allGrSignatures_17.10.2025$minusClusterD_BloodCellsDown

flat_allGrSignatures_17.10.2025$global_GR_genes_globalUp5TissuesDerivedCells


intersect(flat_allGrSignatures_17.10.2025$LungCellsUp,
          flat_allGrSignatures_17.10.2025$LungCellsDown)


intersect(flat_allGrSignatures_17.10.2025$BloodCellsUp,
          flat_allGrSignatures_17.10.2025$BloodCellsDown)

tmp$processed$original_data$df %>% 
  filter(Var1 %in% c("global_GR_genes_globalUp5TissuesDerivedCells",
                     "minusGlobalUp5TissuesDerivedCells_NeuralCellsUp",
                     "minusGlobalUp5TissuesDerivedCells_BloodCellsUp",
                     "minusGlobalUp5TissuesDerivedCells_LungCellsUp")) %>% 
  filter(grepl("Up5", Var2)) %>% 
  filter(Var1 != "global_GR_genes_globalUp5TissuesDerivedCells") %>% 
  filter(Var2 != "global_GR_genes_globalUp5TissuesDerivedCells") %>% 
  filter(!is.na(odds_ratio)) %>% 
  select(-c(fdr_value, fdr)) %>% 
  mutate(pair = pmap_chr(list(Var1, Var2), ~paste(sort(c(...)), collapse = "_vs_"))) %>%
  distinct(pair, .keep_all = TRUE) %>%
  select(-pair)


tmp$processed$original_data$df %>% 
  filter(Var1 %in% c("global_GR_genes_globalDown5TissuesDerivedCells",
                     "minusGlobalDown5TissuesDerivedCells_NeuralCellsDown",
                     "minusGlobalDown5TissuesDerivedCells_BloodCellsDown",
                     "minusGlobalDown5TissuesDerivedCells_LungCellsDown")) %>% 
  filter(grepl("Down5", Var2)) %>% 
  filter(Var1 != "global_GR_genes_globalDown5TissuesDerivedCells") %>% 
  filter(Var2 != "global_GR_genes_globalDown5TissuesDerivedCells") %>% 
  filter(!is.na(odds_ratio)) %>% 
  select(-c(fdr_value, fdr)) %>% 
  mutate(pair = pmap_chr(list(Var1, Var2), ~paste(sort(c(...)), collapse = "_vs_"))) %>%
  distinct(pair, .keep_all = TRUE) %>%
  select(-pair)


tmp$processed$original_data$df %>% 
  filter(Var1 %in% c("global_GR_genes_globalUp4TissuesDerivedCells",
                     "minusGlobalUp4TissuesDerivedCells_NeuralCellsUp",
                     "minusGlobalUp4TissuesDerivedCells_BloodCellsUp",
                     "minusGlobalUp4TissuesDerivedCells_LungCellsUp")) %>% 
  filter(grepl("Up4", Var2)) %>% 
  filter(Var1 != "global_GR_genes_globalUp4TissuesDerivedCells") %>% 
  filter(Var2 != "global_GR_genes_globalUp4TissuesDerivedCells") %>% 
  filter(!is.na(odds_ratio)) %>% 
  select(-c(fdr_value, fdr)) %>% 
  mutate(pair = pmap_chr(list(Var1, Var2), ~paste(sort(c(...)), collapse = "_vs_"))) %>%
  distinct(pair, .keep_all = TRUE) %>%
  select(-pair)


tmp$processed$original_data$df %>% 
  filter(Var1 %in% c("global_GR_genes_globalDown4TissuesDerivedCells",
                     "minusGlobalDown4TissuesDerivedCells_NeuralCellsDown",
                     "minusGlobalDown4TissuesDerivedCells_BloodCellsDown",
                     "minusGlobalDown4TissuesDerivedCells_LungCellsDown")) %>% 
  filter(grepl("Down4", Var2)) %>% 
  filter(Var1 != "global_GR_genes_globalDown4TissuesDerivedCells") %>% 
  filter(Var2 != "global_GR_genes_globalDown4TissuesDerivedCells") %>% 
  filter(!is.na(odds_ratio)) %>% 
  select(-c(fdr_value, fdr)) %>% 
  mutate(pair = pmap_chr(list(Var1, Var2), ~paste(sort(c(...)), collapse = "_vs_"))) %>%
  distinct(pair, .keep_all = TRUE) %>%
  select(-pair)





tmp$processed$original_data$df %>% 
  filter(Var1 %in% c("global_GR_genes_globalUp6TissuesDerivedCells",
                     "minusGlobalUp6TissuesDerivedCells_NeuralCellsUp",
                     "minusGlobalUp6TissuesDerivedCells_BloodCellsUp",
                     "minusGlobalUp6TissuesDerivedCells_LungCellsUp")) %>% 
  filter(grepl("Up6", Var2)) %>% 
  filter(Var1 != "global_GR_genes_globalUp6TissuesDerivedCells") %>% 
  filter(Var2 != "global_GR_genes_globalUp6TissuesDerivedCells") %>% 
  filter(!is.na(odds_ratio)) %>% 
  select(-c(fdr_value, fdr)) %>% 
  mutate(pair = pmap_chr(list(Var1, Var2), ~paste(sort(c(...)), collapse = "_vs_"))) %>%
  distinct(pair, .keep_all = TRUE) %>%
  select(-pair)


tmp$processed$original_data$df %>% 
  filter(Var1 %in% c("global_GR_genes_globalDown6TissuesDerivedCells",
                     "minusGlobalDown6TissuesDerivedCells_NeuralCellsDown",
                     "minusGlobalDown6TissuesDerivedCells_BloodCellsDown",
                     "minusGlobalDown6TissuesDerivedCells_LungCellsDown")) %>% 
  filter(grepl("Down6", Var2)) %>% 
  filter(Var1 != "global_GR_genes_globalDown6TissuesDerivedCells") %>% 
  filter(Var2 != "global_GR_genes_globalDown6TissuesDerivedCells") %>% 
  filter(!is.na(odds_ratio)) %>% 
  select(-c(fdr_value, fdr)) %>% 
  mutate(pair = pmap_chr(list(Var1, Var2), ~paste(sort(c(...)), collapse = "_vs_"))) %>%
  distinct(pair, .keep_all = TRUE) %>%
  select(-pair)


tmp$processed$original_data$df %>% 
  filter(Var1 %in% c("minusClustersKPO_NeuralCellsUp",
                     "minusClustersKPO_BloodCellsUp",
                     "minusClustersKPO_LungCellsUp")) %>% 
  filter(grepl("minusClustersKPO_", Var2)) %>% 
  filter(!is.na(odds_ratio)) %>% 
  select(-c(fdr_value, fdr)) %>% 
  mutate(pair = pmap_chr(list(Var1, Var2), ~paste(sort(c(...)), collapse = "_vs_"))) %>%
  distinct(pair, .keep_all = TRUE) %>%
  select(-pair)

tmp$processed$original_data$df %>% 
  filter(Var1 %in% c("minusClusterD_NeuralCellsDown",
                     "minusClusterD_BloodCellsDown",
                     "minusClusterD_LungCellsDown")) %>% 
  filter(grepl("minusClusterD_", Var2)) %>% 
  filter(!is.na(odds_ratio)) %>% 
  select(-c(fdr_value, fdr)) %>% 
  mutate(pair = pmap_chr(list(Var1, Var2), ~paste(sort(c(...)), collapse = "_vs_"))) %>%
  distinct(pair, .keep_all = TRUE) %>%
  select(-pair)


# ##############################################################################

tmp$processed$original_data$df %>% 
  filter(Var1 %in% c("NeuralCellsUp",
                     "BloodCellsUp",
                     "LungCellsUp",
                     "NeuralCellsDown",
                     "BloodCellsDown",
                     "LungCellsDown",
                     "global_GR_genes_globalUp5TissuesDerivedCells",
                     "global_GR_genes_globalDown5TissuesDerivedCells")) %>% 
  filter(Var2 %in% c("NeuralCellsUp",
                     "BloodCellsUp",
                     "LungCellsUp",
                     "NeuralCellsDown",
                     "BloodCellsDown",
                     "LungCellsDown",
                     "global_GR_genes_globalUp5TissuesDerivedCells",
                     "global_GR_genes_globalDown5TissuesDerivedCells")) %>% 
  filter(!is.na(odds_ratio)) %>% 
  select(-c(fdr_value, fdr)) %>% 
  mutate(pair = pmap_chr(list(Var1, Var2), ~paste(sort(c(...)), collapse = "_vs_"))) %>%
  distinct(pair, .keep_all = TRUE) %>%
  select(-pair)


