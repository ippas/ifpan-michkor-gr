
gene_list <- pgc_annotation_geneCenter50kb_p1e4 %>% 
  filter(pvalue < 0.0001) %>% 
  
  select(gene_symbol, source_file) %>%
  distinct() %>%
  group_by(source_file) %>%
  summarise(genes = list(unique(gene_symbol))) %>%
  deframe()

run_full_overlap_analysis(
  gene_lists = c(flat_allGrSignatures_17.10.2025[c("minusGlobalUp5TissuesDerivedCells_BloodCellsUp",
                                                   "minusGlobalUp5TissuesDerivedCells_LungCellsUp",
                                                   "minusGlobalUp5TissuesDerivedCells_NeuralCellsUp",
                                                   "minusGlobalDown5TissuesDerivedCells_BloodCellsDown",
                                                   "minusGlobalDown5TissuesDerivedCells_LungCellsDown",
                                                   "minusGlobalDown5TissuesDerivedCells_NeuralCellsDown",
                                                   "global_GR_genes_globalUp5TissuesDerivedCells",
                                                   "global_GR_genes_globalDown5TissuesDerivedCells"
  )], gene_list),
  total_genes = hgnc_symbols_vector_v110,
  cols_to_filter = c("minusGlobalUp5TissuesDerivedCells_BloodCellsUp",
    "minusGlobalUp5TissuesDerivedCells_LungCellsUp",
    "minusGlobalUp5TissuesDerivedCells_NeuralCellsUp",
    "minusGlobalDown5TissuesDerivedCells_BloodCellsDown",
    "minusGlobalDown5TissuesDerivedCells_LungCellsDown",
    "minusGlobalDown5TissuesDerivedCells_NeuralCellsDown",
    "global_GR_genes_globalUp5TissuesDerivedCells",
    "global_GR_genes_globalDown5TissuesDerivedCells"
  ),
  rows_to_filter = names(gene_list),
  plot_title_or = "",
  triangle_mode = "full",
  fdr_threshold = 1,
  data_type = "original_data",
  verbose = F
) -> pgcGrSignatures_overlapChi2


# pgcGrSignatures_overlapChi2$processed$significant_data$df %>% dim
# 
# pgcGrSignatures_overlapChi2$processed$original_data$df %>% dim

tmp$processed$original_data$df %>% 
  filter(
    Var1 %in% names(gene_list),
    Var2 %in% names(flat_allGrSignatures_17.10.2025)
  ) %>%
  select(-fdr) %>% 
  group_by(Var1) %>%
  nest() %>% 
  mutate(
    data = map(data, ~ .x %>% 
                 mutate(fdr = p.adjust(p_value, method = "fdr"))
    )
  ) %>%
  unnest(data) %>% filter(p_value < 0.05) %>% 
  # as.data.frame() %>% filter(fdr < 0.1)
  
  filter(gene_overlap_count > 2) %>% 
  filter(log2_odds_ratio > 0) %>% 
  # filter(p_value < 0.05) %>%
  filter(grepl("5", Var1)) %>%
  # filter(grepl("Blood|Neural", Var1)) %>%
  # filter(!grepl("Cluster", Var1)) %>%
  filter(Var1 %in% c("minusGlobalUp5TissuesDerivedCells_BloodCellsUp",
                     "minusGlobalUp5TissuesDerivedCells_LungCellsUp",
                     "minusGlobalUp5TissuesDerivedCells_NeuralCellsUp",
                     "minusGlobalDown5TissuesDerivedCells_BloodCellsDown",
                     "minusGlobalDown5TissuesDerivedCells_LungCellsDown",
                     "minusGlobalDown5TissuesDerivedCells_NeuralCellsDown",
                     "global_GR_genes_globalUp5TissuesDerivedCells",
                     "global_GR_genes_globalDown5TissuesDerivedCells"
                     )) %>%
  # filter(fdr_value < 0.1) %>% 
  .$overlap_genes %>% strsplit(",") %>% unlist %>% unique() %>% length()



tmp$processed$original_data$df %>% 
  filter(Var1 == "global_GR_genes_globalUp5TissuesDerivedCells" & Var2 == "global_GR_genes_globalDown5TissuesDerivedCells")
