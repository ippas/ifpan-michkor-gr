pgcGrSignatures_overlapChi2$processed$original_data$df %>% 
  filter(Var2 %in% c("minusGlobalUp5TissuesDerivedCells_BloodCellsUp",
                     "minusGlobalUp5TissuesDerivedCells_LungCellsUp",
                     "minusGlobalUp5TissuesDerivedCells_NeuralCellsUp",
                     "minusGlobalDown5TissuesDerivedCells_BloodCellsDown",
                     "minusGlobalDown5TissuesDerivedCells_LungCellsDown",
                     "minusGlobalDown5TissuesDerivedCells_NeuralCellsDown",
                     "global_GR_genes_globalUp5TissuesDerivedCells",
                     "global_GR_genes_globalDown5TissuesDerivedCells"
  )) %>% 
  group_by(Var2) %>%
  nest() %>% 
  mutate(
    data = map(data, ~ .x %>% 
                 mutate(fdr = p.adjust(p_value, method = "fdr"))
    )
  ) %>% 
  unnest(data) %>% 
  filter(gene_overlap_count > 2) %>% 
  filter(log2_odds_ratio > 0) %>%
  filter(p_value < 0.05)
  

read.delim(gzfile("data/PGC/mergedPGC_annotation_geneCenter50kb_p1e4.tsv"),
           header = TRUE, sep = "\t") %>% 
  filter(source_file != "") 

gene_list %>% lapply(., length) %>% unname() %>% unlist %>% summary

pgcGrSignatures_overlapChi2$processed$original_data$df %>% 
  filter(Var2 %in% c("minusGlobalUp5TissuesDerivedCells_BloodCellsUp",
                     "minusGlobalUp5TissuesDerivedCells_LungCellsUp",
                     "minusGlobalUp5TissuesDerivedCells_NeuralCellsUp",
                     "minusGlobalDown5TissuesDerivedCells_BloodCellsDown",
                     "minusGlobalDown5TissuesDerivedCells_LungCellsDown",
                     "minusGlobalDown5TissuesDerivedCells_NeuralCellsDown",
                     "global_GR_genes_globalUp5TissuesDerivedCells",
                     "global_GR_genes_globalDown5TissuesDerivedCells"
  )) %>% 
  group_by(Var2) %>%
  nest() %>% 
  mutate(
    data = map(data, ~ .x %>% 
                 mutate(fdr = p.adjust(p_value, method = "fdr"))
    )
  ) %>% 
  unnest(data) %>% 
  filter(gene_overlap_count > 2) %>% 
  filter(log2_odds_ratio > 0) %>%
  filter(p_value < 0.05) %>% .$Var2 -> pgc_grSignatures_vector_p0.05


gene_list[pgc_phenotypes_vector_p0.05] %>% 
  lapply(., length)







heatmap_overlap_log2OR_ggplot(
  data_list = pgcGrSignatures_overlapChi2$processed,
  data_type = "original_data",
  triangle_mode = "full",
  title = "Odds ratio overlap",
  color_scale_range = c(-3, 3),
  text_contrast_range = c(-30, 4.9),
  palette = c("#c6d3e3", "white", "darkred"),
  rows_to_filter = pgc_phenotypes_vector_p0.05,
  cols_to_filter = pgc_grSignatures_vector_p0.05,
  color_rects = c("#97C426", "#2F4603")
) -> p1


custom_legend <- create_customRect_patch_legend(
  colors = c("#97C426", "#2F4603"),
  labels = c("p < 0.05", "p < 0.01"),
  spacing = 4,
  box_linewidth = 4
)

p1 + custom_legend

wrap_plots(
  p1,
  custom_legend,
  ncol = 1,
  heights = c(10, 1)
)



# prepare table
pgcGrSignatures_overlapChi2$processed$original_data$df %>% 
  filter(Var2 %in% c("minusGlobalUp5TissuesDerivedCells_BloodCellsUp",
                     "minusGlobalUp5TissuesDerivedCells_LungCellsUp",
                     "minusGlobalUp5TissuesDerivedCells_NeuralCellsUp",
                     "minusGlobalDown5TissuesDerivedCells_BloodCellsDown",
                     "minusGlobalDown5TissuesDerivedCells_LungCellsDown",
                     "minusGlobalDown5TissuesDerivedCells_NeuralCellsDown",
                     "global_GR_genes_globalUp5TissuesDerivedCells",
                     "global_GR_genes_globalDown5TissuesDerivedCells"
  )) %>% 
  group_by(Var2) %>%
  nest() %>% 
  mutate(
    data = map(data, ~ .x %>% 
                 mutate(fdr = p.adjust(p_value, method = "fdr"))
    )
  ) %>% 
  unnest(data) %>% 
  group_by(Var2) %>% 
  nest %>% 
  mutate(
    n_genes_all = map(data, ~ .x %>% 
                    filter(overlap_genes != "") %>% 
                    pull(overlap_genes) %>% 
                    strsplit(",") %>% 
                    unlist() %>% 
                    unique() %>% 
                    length())
  ) %>% 
  mutate(
    n_genes_p0.05 = map(data, ~ .x %>% 
                        filter(p_value < 0.05) %>% 
                        filter(log2_odds_ratio > 0) %>%   
                        filter(overlap_genes != "") %>% 
                        filter(gene_overlap_count > 2) %>% 
                        pull(overlap_genes) %>% 
                        strsplit(",") %>% 
                        unlist() %>% 
                        unique() %>% 
                        length())
  ) %>% 
  mutate(
    n_genes_p0.01 = map(data, ~ .x %>% 
                          filter(p_value < 0.01) %>% 
                          filter(log2_odds_ratio > 0) %>%   
                          filter(overlap_genes != "") %>% 
                          filter(gene_overlap_count > 2) %>% 
                          pull(overlap_genes) %>% 
                          strsplit(",") %>% 
                          unlist() %>% 
                          unique() %>% 
                          length())
  ) %>% 
  unnest(n_genes_all, n_genes_p0.05, n_genes_p0.01)
  # filter(gene_overlap_count > 2) %>% 
  # filter(log2_odds_ratio > 0) %>%
  # filter(p_value < 0.05) %>% 
