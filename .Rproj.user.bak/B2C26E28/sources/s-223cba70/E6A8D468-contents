
# top20 brain_up
selected_ranks <- grSignature_genebass_associations_FDRMonteCarlo$original_genebass_association$brain_up %>% 
  mutate(fdr_0.05 = ifelse(pvalue < p5_pvalue, T, F)) %>%
  mutate(fdr_0.1 = ifelse(pvalue < p10_pvalue, T, F)) %>% 
  mutate(fdr_0.2 = ifelse(pvalue < p20_pvalue, T, F)) %>% 
  filter(fdr_0.2) %>% 
  head(20) %>% .$rank

svg("data/genebass/figures/brainUp_genebassSKAT_association_top20_fdr0.2.svg", height = 12, width = 26)

p1 <- plot_permutation_pvalues_by_rank_v3(
  signature_data = grSignature_genebass_associations_FDRMonteCarlo$original_genebass_association$brain_up,
  ranks_to_plot = selected_ranks,
  x_axis_limit = 16,
  pvalue_types_to_plot = c("pvalue", "median_pvalue", "p5_pvalue", "p10_pvalue", "p20_pvalue", "min_pvalue", "max_pvalue"),
  title_text_size = 20,
  subtitle_text_size = 16,
  axis_title_size = 18,
  axis_text_x_size = 18,
  axis_text_y_size = 18,
  legend_text_size = 18,
  legend_title_size = 18,
  use_mono_font = FALSE,
  point_size = 3,
  pvalue_color_map_override = c(
    "pvalue" = "blue",
    "median_pvalue" = "#f0c1c1",
    "p5_pvalue" = "#4b0b0b",
    "p10_pvalue" = "#6b1e1e",
    "p20_pvalue" = "#9c5757",
    "min_pvalue" = "#558b2f",
    "max_pvalue" = "#aed581"
  )
)

dev.off()

# top20 metasignature_down
selected_ranks <- grSignature_genebass_associations_FDRMonteCarlo$original_genebass_association$metasignature_down %>% 
  mutate(fdr_0.05 = ifelse(pvalue < p5_pvalue, T, F)) %>%
  mutate(fdr_0.1 = ifelse(pvalue < p10_pvalue, T, F)) %>% 
  mutate(fdr_0.2 = ifelse(pvalue < p20_pvalue, T, F)) %>% 
  filter(fdr_0.2) %>% 
  head(20) %>% .$rank


svg("data/genebass/figures/metasignatureDown_genebassSKAT_association_top20_fdr0.2.svg", height = 12, width = 26)
p2 <- plot_permutation_pvalues_by_rank_v3(
  signature_data = grSignature_genebass_associations_FDRMonteCarlo$original_genebass_association$metasignature_down,
  ranks_to_plot = selected_ranks,
  x_axis_limit = 16,
  pvalue_types_to_plot = c("pvalue", "median_pvalue", "p5_pvalue", "p10_pvalue", "p20_pvalue", "min_pvalue", "max_pvalue"),
  title_text_size = 20,
  subtitle_text_size = 16,
  axis_title_size = 18,
  axis_text_x_size = 18,
  axis_text_y_size = 18,
  legend_text_size = 18,
  legend_title_size = 18,
  use_mono_font = FALSE,
  point_size = 3,
  pvalue_color_map_override = c(
    "pvalue" = "blue",
    "median_pvalue" = "#f0c1c1",
    "p5_pvalue" = "#4b0b0b",
    "p10_pvalue" = "#6b1e1e",
    "p20_pvalue" = "#9c5757",
    "min_pvalue" = "#558b2f",
    "max_pvalue" = "#aed581"
  )
)
dev.off()

selected_ranks <- grSignature_genebass_associations_FDRMonteCarlo$original_genebass_association$metasignature_up %>% 
  mutate(fdr_0.05 = ifelse(pvalue < p5_pvalue, T, F)) %>%
  mutate(fdr_0.1 = ifelse(pvalue < p10_pvalue, T, F)) %>% 
  mutate(fdr_0.2 = ifelse(pvalue < p20_pvalue, T, F)) %>% 
  filter(fdr_0.2) %>% 
  head(20) %>% .$rank


svg("data/genebass/figures/metasignatureUp_genebassSKAT_association_top20_fdr0.2.svg", height = 12, width = 26)
p3 <- plot_permutation_pvalues_by_rank_v3(
  signature_data = grSignature_genebass_associations_FDRMonteCarlo$original_genebass_association$metasignature_up,
  ranks_to_plot = selected_ranks,
  x_axis_limit = 16,
  pvalue_types_to_plot = c("pvalue", "median_pvalue", "p5_pvalue", "p10_pvalue", "p20_pvalue", "min_pvalue", "max_pvalue"),
  title_text_size = 20,
  subtitle_text_size = 16,
  axis_title_size = 18,
  axis_text_x_size = 18,
  axis_text_y_size = 18,
  legend_text_size = 18,
  legend_title_size = 18,
  use_mono_font = FALSE,
  point_size = 3,
  pvalue_color_map_override = c(
    "pvalue" = "blue",
    "median_pvalue" = "#f0c1c1",
    "p5_pvalue" = "#4b0b0b",
    "p10_pvalue" = "#6b1e1e",
    "p20_pvalue" = "#9c5757",
    "min_pvalue" = "#558b2f",
    "max_pvalue" = "#aed581"
  )
)
dev.off()


p1 / p2 / p3

svg("data/genebass/figures/GRsignaturesCombine_genebassSKAT_association_top20_fdr0.2.svg", height = 36, width = 28)
p3 / p2 / p1
dev.off()
