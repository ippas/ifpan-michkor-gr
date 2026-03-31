

plot_icd10_bar(
  preprocessing_summaryTable_grSystemicTissues_icd10F_df,
  fill_colors = c(
    lung     = "#E09F3E",  
    blood    = "#9E2A2B",
    neural   = "#540B0E",
    systemic = "#335C67"),   # przygaszona czerwień),
  mean_signif_cs,
  border_color = "black",
  # x_range = c(0,5),
  y_order = rev(c("F3x", "F2x", "F0x", "F1x", "F4x", "F8x", "F5x", "F6x", "F7x", "F9x"))
) -> p1

plot_icd10_bar(
  preprocessing_summaryTable_grSystemicTissues_icd10F_df,
  fill_colors = c(
    lung     = "#540B0E",  
    blood    = "#9E2A2B",
    neural   = "#E09F3E",
    systemic = "#335C67"),   # przygaszona czerwień),
  mean_signif_cs,
  border_color = "black",
  # x_range = c(0,5),
  y_order = rev(c("F3x", "F2x", "F0x", "F1x", "F4x", "F8x", "F5x", "F6x", "F7x", "F9x"))
) -> p2


plot_icd10_bar(
  preprocessing_summaryTable_grSystemicTissues_icd10F_df,
  fill_colors = c(
    lung     = "#540B0E",  
    blood    = "#9E2A2B",
    neural   = "#FFF3B0",
    systemic = "#335C67"),   # przygaszona czerwień),
  mean_signif_cs,
  border_color = "black",
  # x_range = c(0,5),
  y_order = rev(c("F3x", "F2x", "F0x", "F1x", "F4x", "F8x", "F5x", "F6x", "F7x", "F9x"))
) -> p3

plot_icd10_bar(
  preprocessing_summaryTable_grSystemicTissues_icd10F_df,
  fill_colors = c(
    lung     = "#9E2A2B",  
    blood    = "#E09F3E",
    neural   = "#FFF3B0",
    systemic = "#335C67"),   # przygaszona czerwień),
  mean_signif_cs,
  border_color = "black",
  # x_range = c(0,5),
  y_order = rev(c("F3x", "F2x", "F0x", "F1x", "F4x", "F8x", "F5x", "F6x", "F7x", "F9x"))
) -> p4

plot_icd10_bar(
  preprocessing_summaryTable_grSystemicTissues_icd10F_df,
  fill_colors = c(
    lung     = "#FFF3B0",  
    blood    = "#E09F3E",
    neural   = "#9E2A2B",
    systemic = "#335C67"),   # przygaszona czerwień),
  mean_signif_cs,
  border_color = "black",
  # x_range = c(0,5),
  y_order = rev(c("F3x", "F2x", "F0x", "F1x", "F4x", "F8x", "F5x", "F6x", "F7x", "F9x"))
) -> p5





p2 + p4  + p1 + p3 +  p5 




plot_icd10_bar(
  preprocessing_summaryTable_grSystemicTissues_icd10F_df,
  fill_colors = c(
    lung     = "#540B0E",  
    blood    = "#9E2A2B",
    neural   = "#E09F3E",
    systemic = "#335C67"),   # przygaszona czerwień),
  mean_signif_cs,
  border_color = "black",
  # x_range = c(0,5),
  y_order = rev(c("F3x", "F2x", "F0x", "F1x", "F4x", "F8x", "F5x", "F6x", "F7x", "F9x"))
) -> p1


plot_icd10_bar(
  preprocessing_summaryTable_grSystemicTissues_icd10F_df %>% 
    mutate(mean_signif_cs = ifelse(mean_signif_cs == 0, mean_signif_cs + 0.008, mean_signif_cs)),
  fill_colors = c(
    lung     = "#540B0E",  
    blood    = "#9E2A2B",
    neural   = "#E09F3E",
    systemic = "#335C67"),   # przygaszona czerwień),
  mean_signif_cs,
  # border_color = "black",
  # x_range = c(0,5),
  y_order = rev(c("F3x", "F2x", "F0x", "F1x", "F4x", "F8x", "F5x", "F6x", "F7x", "F9x"))
) -> p2

p1 + p2


svg("/home/mateusz/projects/ifpan-michkor-gr/results_v2/overlap/overlap_grTissues_GWASCatalogDisGeNETgenebass/figures/barplot_grTissuesSystemic_perTissue_meanCS_v3_05.03.2026.svg",
    width = 8,
    height = 8)

(p1 + p2 + p1 + p2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

dev.off()

