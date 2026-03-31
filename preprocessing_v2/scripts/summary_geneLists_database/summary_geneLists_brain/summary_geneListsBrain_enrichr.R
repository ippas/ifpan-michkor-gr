
# ##############################################################################
# ---- filtered ----
# ##############################################################################
plot_integrated_enrichment_stacked(
  gene_summary_subset = gene_summary_list$n10_short_time,
  subdata = "summary_up",
  enrichr_source = "enrichr_3pub",
  fdr_threshold = 0.05,
  y_limits = c(0, 45),
  axis_text_size = 24,
  axis_title_size = 28,
  label_text_size = 22,
  legend_text_size = 28,
  x_text_angle = 25
)

plot_integrated_enrichment_stacked(
  gene_summary_subset = gene_summary_list$n10_short_time,
  subdata = "summary_down",
  enrichr_source = "enrichr_3pub",
  fdr_threshold = 0.05,
  y_limits = c(0, 8),
  palette_sources = c(
    "ChEA_2022" = "#08306b",
    "LINCS_L1000" = "#2171b5",
    "CellMarker_2024" = "#bdd7e7"
  ),
  axis_text_size = 24,
  axis_title_size = 28,
  label_text_size = 22,
  legend_text_size = 28,
  x_text_angle = 25
)





plot_integrated_enrichment_stacked(
  gene_summary_subset = gene_summary_list$n10_short_time,
  subdata = "summary_up",
  enrichr_source = "enrichr_4pub",
  fdr_threshold = 0.05,
  y_limits = c(0, 25),
  axis_text_size = 24,
  axis_title_size = 28,
  label_text_size = 22,
  legend_text_size = 28,
  x_text_angle = 25
)

plot_integrated_enrichment_stacked(
  gene_summary_subset = gene_summary_list$n10_short_time,
  subdata = "summary_down",
  enrichr_source = "enrichr_4pub",
  fdr_threshold = 0.05,
  y_limits = c(0, 3),
  palette_sources = c(
    "ChEA_2022" = "#08306b",
    "LINCS_L1000" = "#2171b5",
    "CellMarker_2024" = "#bdd7e7"
  ),
  axis_text_size = 24,
  axis_title_size = 28,
  label_text_size = 22,
  legend_text_size = 28,
  x_text_angle = 25
)




plot_integrated_enrichment_stacked(
  gene_summary_subset = gene_summary_list$n10_short_time,
  subdata = "summary_up",
  enrichr_source = "enrichr_5pub",
  fdr_threshold = 0.05,
  y_limits = c(0, 8),
  axis_text_size = 24,
  axis_title_size = 28,
  label_text_size = 22,
  legend_text_size = 28,
  x_text_angle = 25
)

plot_integrated_enrichment_stacked(
  gene_summary_subset = gene_summary_list$n10_short_time,
  subdata = "summary_down",
  enrichr_source = "enrichr_5pub",
  fdr_threshold = 0.05,
  y_limits = c(0, 3),
  palette_sources = c(
    "ChEA_2022" = "#08306b",
    "LINCS_L1000" = "#2171b5",
    "CellMarker_2024" = "#bdd7e7"
  ),
  axis_text_size = 24,
  axis_title_size = 28,
  label_text_size = 22,
  legend_text_size = 28,
  x_text_angle = 25
)






# ##############################################################################
# ---- filtered ----
# ##############################################################################

plot_integrated_enrichment_stacked(
  gene_summary_subset = gene_summary_list$all,
  subdata = "summary_up",
  enrichr_source = "enrichr_3pub",
  fdr_threshold = 0.05,
  y_limits = c(0, 45),
  axis_text_size = 24,
  axis_title_size = 28,
  label_text_size = 22,
  legend_text_size = 28,
  x_text_angle = 25
)

plot_integrated_enrichment_stacked(
  gene_summary_subset = gene_summary_list$all,
  subdata = "summary_down",
  enrichr_source = "enrichr_3pub",
  fdr_threshold = 0.05,
  y_limits = c(0, 10),
  palette_sources = c(
    "ChEA_2022" = "#08306b",
    "LINCS_L1000" = "#2171b5",
    "CellMarker_2024" = "#bdd7e7"
  ),
  axis_text_size = 24,
  axis_title_size = 28,
  label_text_size = 22,
  legend_text_size = 28,
  x_text_angle = 25
)



plot_integrated_enrichment_stacked(
  gene_summary_subset = gene_summary_list$all,
  subdata = "summary_up",
  enrichr_source = "enrichr_4pub",
  fdr_threshold = 0.05,
  y_limits = c(0, 30),
  axis_text_size = 24,
  axis_title_size = 28,
  label_text_size = 22,
  legend_text_size = 28,
  x_text_angle = 25
)

plot_integrated_enrichment_stacked(
  gene_summary_subset = gene_summary_list$all,
  subdata = "summary_down",
  enrichr_source = "enrichr_4pub",
  fdr_threshold = 0.05,
  y_limits = c(0, 4),
  palette_sources = c(
    "ChEA_2022" = "#08306b",
    "LINCS_L1000" = "#2171b5",
    "CellMarker_2024" = "#bdd7e7"
  ),
  axis_text_size = 24,
  axis_title_size = 28,
  label_text_size = 22,
  legend_text_size = 28,
  x_text_angle = 25
)


plot_integrated_enrichment_stacked(
  gene_summary_subset = gene_summary_list$all,
  subdata = "summary_up",
  enrichr_source = "enrichr_5pub",
  fdr_threshold = 0.05,
  y_limits = c(0, 12),
  axis_text_size = 24,
  axis_title_size = 28,
  label_text_size = 22,
  legend_text_size = 28,
  x_text_angle = 25
)

plot_integrated_enrichment_stacked(
  gene_summary_subset = gene_summary_list$all,
  subdata = "summary_down",
  enrichr_source = "enrichr_5pub",
  fdr_threshold = 0.05,
  y_limits = c(0, 2),
  palette_sources = c(
    "ChEA_2022" = "#08306b",
    "LINCS_L1000" = "#2171b5",
    "CellMarker_2024" = "#bdd7e7"
  ),
  axis_text_size = 24,
  axis_title_size = 28,
  label_text_size = 22,
  legend_text_size = 28,
  x_text_angle = 25
)



gene_summary_list$n10_short_time$summary_up$enrichr_summary_source$enrichr_3pub$CellMarker_2024 %>% head(10) %>% 
  select(-c(Old.P.value, Old.Adjusted.P.value))

gene_summary_list$n10_short_time$summary_up$enrichr_summary_source$enrichr_3pub$ChEA_2022 %>% head(10) %>% 
  select(-c(Old.P.value, Old.Adjusted.P.value))

gene_summary_list$n10_short_time$summary_up$enrichr_summary_source$enrichr_3pub$LINCS_L1000_Chem_Pert_Consensus_Sigs %>% head(10) %>% .[, 1:8] %>% 
  select(-c(Old.P.value, Old.Adjusted.P.value))


gene_summary_list$all$summary_up$enrichr_summary_source$enrichr_3pub$CellMarker_2024 %>% head


intersect(gene_summary_list$all$summary_up$enrichr_summary_source$genes_3pub,
          gene_summary_list$all$summary_down$enrichr_summary_source$genes_3pub
          )
gene_summary_list$all$summary_up$enrichr_summary_source$genes_3pub

gene_summary_list$all$summary_down$enrichr_summary_source$genes_3pub


