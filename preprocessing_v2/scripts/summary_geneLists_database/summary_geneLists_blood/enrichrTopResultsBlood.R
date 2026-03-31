# ##############################################################################
# ---- uses data ----
# ##############################################################################

gene_summaryList_blood




# ##############################################################################


gene_summaryList_blood$n10_short_time$summary_up$enrichr_summary_source$enrichr_3pub$CellMarker_2024 %>% 
  plot_enrichr_results_logp_clean_tf 



# "ChEA_2022" = "#4b0b0b",
# "LINCS_L1000" = "#9c5757",
# "CellMarker_2024" = "#f0c1c1"

# "ChEA_2022" = "#08306b",
# "LINCS_L1000" = "#2171b5",
# "CellMarker_2024" = "#bdd7e7"

# ##############################################################################
# ---- CellMarker ----
# ##############################################################################

plot_enrichr_barplot_topResults(
  data = gene_summaryList_blood$n10_short_time$summary_up$enrichr_summary_source$enrichr_3pub$CellMarker_2024,
  x_axis = "Adjusted.P.value",
  y_axis = "Term",
  y_label = "Top Enriched Cell Types (blood - up)",
  x_label = expression(-log[10](FDR)),
  top_n = 20,
  label_hjust = -0.1,
  axis_label_size = 18,
  axis_text_size = 14,
  label_size = 12,
  fill_color =  "#f0c1c1"
) -> p1

plot_enrichr_barplot_topResults(
  data = gene_summaryList_blood$n10_short_time$summary_down$enrichr_summary_source$enrichr_3pub$CellMarker_2024,
  x_axis = "Adjusted.P.value",
  y_axis = "Term",
  y_label = "Top Enriched Cell Types (blood - down)",
  top_n = 20,
  label_hjust = -0.1,
  axis_label_size = 18,
  axis_text_size = 14,
  label_size = 12,
  fill_color =  "#bdd7e7",
  x_label = expression(-log[10](FDR)),
) -> p2

p1 / p2


plot_enrichr_barplot_topResults(
  data = gene_summaryList_blood$n10_short_time$summary_up$enrichr_summary_source$enrichr_4pub$CellMarker_2024,
  x_axis = "Adjusted.P.value",
  y_axis = "Term",
  y_label = "Top Enriched Cell Types (blood - up)",
  x_label = expression(-log[10](FDR)),
  top_n = 20,
  label_hjust = -0.1,
  axis_label_size = 18,
  axis_text_size = 14,
  label_size = 12,
  fill_color =  "#f0c1c1"
) -> p1

plot_enrichr_barplot_topResults(
  data = gene_summaryList_blood$n10_short_time$summary_down$enrichr_summary_source$enrichr_4pub$CellMarker_2024,
  x_axis = "Adjusted.P.value",
  y_axis = "Term",
  y_label = "Top Enriched Cell Types (blood - down)",
  top_n = 20,
  label_hjust = -0.1,
  axis_label_size = 18,
  axis_text_size = 14,
  label_size = 12,
  fill_color =  "#bdd7e7",
  x_label = expression(-log[10](FDR)),
) -> p2

p1 / p2

plot_enrichr_barplot_topResults(
  data = gene_summaryList_blood$n10_short_time$summary_up$enrichr_summary_source$enrichr_5pub$CellMarker_2024,
  x_axis = "Adjusted.P.value",
  y_axis = "Term",
  y_label = "Top Enriched Cell Types (blood - up)",
  x_label = expression(-log[10](FDR)),
  top_n = 20,
  label_hjust = -0.1,
  axis_label_size = 18,
  axis_text_size = 14,
  label_size = 12,
  fill_color =  "#f0c1c1"
) -> p1

plot_enrichr_barplot_topResults(
  data = gene_summaryList_blood$n10_short_time$summary_down$enrichr_summary_source$enrichr_5pub$CellMarker_2024,
  x_axis = "Adjusted.P.value",
  y_axis = "Term",
  y_label = "Top Enriched Cell Types (blood - down)",
  top_n = 20,
  label_hjust = -0.1,
  axis_label_size = 18,
  axis_text_size = 14,
  label_size = 12,
  fill_color =  "#bdd7e7",
  x_label = expression(-log[10](FDR)),
) -> p2

p1 / p2

# ##############################################################################
# ---- ChEA 2022 ----
# ##############################################################################


plot_enrichr_barplot_topResults(
  data = gene_summaryList_blood$n10_short_time$summary_up$enrichr_summary_source$enrichr_3pub$ChEA_2022,
  x_axis = "Adjusted.P.value",
  y_axis = "Term",
  y_label = "Top Enriched Cell Types (blood - up)",
  x_label = expression(-log[10](FDR)),
  top_n = 20,
  label_hjust = -0.1,
  axis_label_size = 18,
  axis_text_size = 14,
  label_size = 12,
  fill_color =  "#4b0b0b"
) -> p1

plot_enrichr_barplot_topResults(
  data = gene_summaryList_blood$n10_short_time$summary_down$enrichr_summary_source$enrichr_3pub$ChEA_2022,
  x_axis = "Adjusted.P.value",
  y_axis = "Term",
  y_label = "Top Enriched Cell Types (blood - down)",
  top_n = 20,
  label_hjust = -0.1,
  axis_label_size = 18,
  axis_text_size = 14,
  label_size = 12,
  fill_color =  "#08306b",
  x_label = expression(-log[10](FDR)),
) -> p2

p1 / p2


plot_enrichr_barplot_topResults(
  data = gene_summaryList_blood$n10_short_time$summary_up$enrichr_summary_source$enrichr_4pub$ChEA_2022,
  x_axis = "Adjusted.P.value",
  y_axis = "Term",
  y_label = "Top Enriched Cell Types (blood - up)",
  x_label = expression(-log[10](FDR)),
  top_n = 20,
  label_hjust = -0.1,
  axis_label_size = 18,
  axis_text_size = 14,
  label_size = 12,
  fill_color =  "#4b0b0b"
) -> p1

plot_enrichr_barplot_topResults(
  data = gene_summaryList_blood$n10_short_time$summary_down$enrichr_summary_source$enrichr_4pub$ChEA_2022,
  x_axis = "Adjusted.P.value",
  y_axis = "Term",
  y_label = "Top Enriched Cell Types (blood - down)",
  top_n = 20,
  label_hjust = -0.1,
  axis_label_size = 18,
  axis_text_size = 14,
  label_size = 12,
  fill_color =  "#08306b",
  x_label = expression(-log[10](FDR)),
) -> p2

p1 / p2

plot_enrichr_barplot_topResults(
  data = gene_summaryList_blood$n10_short_time$summary_up$enrichr_summary_source$enrichr_5pub$ChEA_2022,
  x_axis = "Adjusted.P.value",
  y_axis = "Term",
  y_label = "Top Enriched Cell Types (blood - up)",
  x_label = expression(-log[10](FDR)),
  top_n = 20,
  label_hjust = -0.1,
  axis_label_size = 18,
  axis_text_size = 14,
  label_size = 12,
  fill_color =  "#4b0b0b"
) -> p1

plot_enrichr_barplot_topResults(
  data = gene_summaryList_blood$n10_short_time$summary_down$enrichr_summary_source$enrichr_5pub$ChEA_2022,
  x_axis = "Adjusted.P.value",
  y_axis = "Term",
  y_label = "Top Enriched Cell Types (blood - down)",
  top_n = 20,
  label_hjust = -0.1,
  axis_label_size = 18,
  axis_text_size = 14,
  label_size = 12,
  fill_color =  "#08306b",
  x_label = expression(-log[10](FDR)),
) -> p2

p1 / p2



# ##############################################################################
# ---- LINCS L1000 ----
# ##############################################################################

plot_enrichr_barplot_topResults(
  data = gene_summaryList_blood$n10_short_time$summary_up$enrichr_summary_source$enrichr_3pub$LINCS_L1000_Chem_Pert_Consensus_Sigs,
  x_axis = "Adjusted.P.value",
  y_axis = "Term",
  y_label = "Top Enriched Cell Types (blood - up)",
  x_label = expression(-log[10](FDR)),
  top_n = 20,
  label_hjust = -0.1,
  axis_label_size = 18,
  axis_text_size = 14,
  label_size = 12,
  fill_color =  "#9c5757"
) -> p1

plot_enrichr_barplot_topResults(
  data = gene_summaryList_blood$n10_short_time$summary_down$enrichr_summary_source$enrichr_3pub$LINCS_L1000_Chem_Pert_Consensus_Sigs,
  x_axis = "Adjusted.P.value",
  y_axis = "Term",
  y_label = "Top Enriched Cell Types (blood - down)",
  top_n = 20,
  label_hjust = -0.1,
  axis_label_size = 18,
  axis_text_size = 14,
  label_size = 12,
  fill_color =  "#2171b5",
  x_label = expression(-log[10](FDR)),
) -> p2

p1 / p2


plot_enrichr_barplot_topResults(
  data = gene_summaryList_blood$n10_short_time$summary_up$enrichr_summary_source$enrichr_4pub$LINCS_L1000_Chem_Pert_Consensus_Sigs,
  x_axis = "Adjusted.P.value",
  y_axis = "Term",
  y_label = "Top Enriched Cell Types (blood - up)",
  x_label = expression(-log[10](FDR)),
  top_n = 20,
  label_hjust = -0.1,
  axis_label_size = 18,
  axis_text_size = 14,
  label_size = 12,
  fill_color =  "#9c5757"
) -> p1

plot_enrichr_barplot_topResults(
  data = gene_summaryList_blood$n10_short_time$summary_down$enrichr_summary_source$enrichr_4pub$LINCS_L1000_Chem_Pert_Consensus_Sigs,
  x_axis = "Adjusted.P.value",
  y_axis = "Term",
  y_label = "Top Enriched Cell Types (blood - down)",
  top_n = 20,
  label_hjust = -0.1,
  axis_label_size = 18,
  axis_text_size = 14,
  label_size = 12,
  fill_color =  "#2171b5",
  x_label = expression(-log[10](FDR)),
) -> p2

p1 / p2


plot_enrichr_barplot_topResults(
  data = gene_summaryList_blood$n10_short_time$summary_up$enrichr_summary_source$enrichr_5pub$LINCS_L1000_Chem_Pert_Consensus_Sigs,
  x_axis = "Adjusted.P.value",
  y_axis = "Term",
  y_label = "Top Enriched Cell Types (blood - up)",
  x_label = expression(-log[10](FDR)),
  top_n = 20,
  label_hjust = -0.1,
  axis_label_size = 18,
  axis_text_size = 14,
  label_size = 12,
  fill_color =  "#9c5757"
) -> p1

plot_enrichr_barplot_topResults(
  data = gene_summaryList_blood$n10_short_time$summary_down$enrichr_summary_source$enrichr_5pub$LINCS_L1000_Chem_Pert_Consensus_Sigs,
  x_axis = "Adjusted.P.value",
  y_axis = "Term",
  y_label = "Top Enriched Cell Types (blood - down)",
  top_n = 20,
  label_hjust = -0.1,
  axis_label_size = 18,
  axis_text_size = 14,
  label_size = 12,
  fill_color =  "#2171b5",
  x_label = expression(-log[10](FDR)),
) -> p2

p1 / p2
