fdr_threshold <- 0.05
gene_set_name = "genes_5pub"       # 👈 nowy argument
enrichr_set_name = "enrichr_5pub"   # 👈 nowy argument


p1 <- plot_AggregateEnrichr_summarySignatures_normalized(
  data = gene_summaryList_brain$n10_short_time,
  subdata = "summary_up",
  x_label = "",
  y_label = "Percent of genes from signature",
  palette_sources = c(
    "ChEA_2022" = "#4b0b0b",
    "LINCS_L1000" = "#9c5757",
    "CellMarker_2024" = "#f0c1c1"
  ),
  axis_text_size = 20,
  axis_title_size = 20,
  legend_text_size = 18,
  signature_text_size = 20,
  fdr_threshold = fdr_threshold,
  number_association = FALSE,
  gene_set_name = gene_set_name,
  enrichr_set_name = enrichr_set_name,
  y_limits = c(0, 100)  # opcjonalnie – proporcje
)

p2 <- plot_AggregateEnrichr_summarySignatures_normalized(
  data = gene_summaryList_brain$n10_short_time,
  subdata = "summary_down",
  # plot_title = "Σ GR-related aggregated enrichments (Blood, normalized, down)",
  x_label = "",
  y_label = "Percent of genes from signature",
  palette_sources = c(
    "ChEA_2022" = "#08306b",
    "LINCS_L1000" = "#2171b5",
    "CellMarker_2024" = "#bdd7e7"
  ),
  axis_text_size = 20,
  axis_title_size = 20,
  legend_text_size = 18,
  signature_text_size = 20,
  fdr_threshold = fdr_threshold,
  number_association = FALSE,
  gene_set_name = gene_set_name,
  enrichr_set_name = enrichr_set_name,
  y_limits = c(0, 100)  # odpowiadające proporcjonalnemu zakresowi
)


p3 <- plot_AggregateEnrichr_summarySignatures_normalized(
  data = gene_summaryList_blood$n10_short_time,
  subdata = "summary_up",
  x_label = "",
  y_label = "Percent of genes from signature",
  palette_sources = c(
    "ChEA_2022" = "#4b0b0b",
    "LINCS_L1000" = "#9c5757",
    "CellMarker_2024" = "#f0c1c1"
  ),
  axis_text_size = 20,
  axis_title_size = 20,
  legend_text_size = 18,
  signature_text_size = 20,
  fdr_threshold = fdr_threshold,
  number_association = FALSE,
  gene_set_name = gene_set_name,
  enrichr_set_name = enrichr_set_name,
  y_limits = c(0, 100)  # opcjonalnie – proporcje
)

p4 <- plot_AggregateEnrichr_summarySignatures_normalized(
  data = gene_summaryList_blood$n10_short_time,
  subdata = "summary_down",
  # plot_title = "Σ GR-related aggregated enrichments (Blood, normalized, down)",
  x_label = "",
  y_label = "Percent of genes from signature",
  palette_sources = c(
    "ChEA_2022" = "#08306b",
    "LINCS_L1000" = "#2171b5",
    "CellMarker_2024" = "#bdd7e7"
  ),
  axis_text_size = 20,
  axis_title_size = 20,
  legend_text_size = 18,
  signature_text_size = 20,
  fdr_threshold = fdr_threshold,
  number_association = FALSE,
  gene_set_name = gene_set_name,
  enrichr_set_name = enrichr_set_name,
  y_limits = c(0, 100)  # odpowiadające proporcjonalnemu zakresowi
)


p5 <- plot_AggregateEnrichr_summarySignatures_normalized(
  data = gene_summaryList_lung$n10_short_time,
  subdata = "summary_up",
  x_label = "",
  y_label = "Percent of genes from signature",
  palette_sources = c(
    "ChEA_2022" = "#4b0b0b",
    "LINCS_L1000" = "#9c5757",
    "CellMarker_2024" = "#f0c1c1"
  ),
  axis_text_size = 20,
  axis_title_size = 20,
  legend_text_size = 18,
  signature_text_size = 20,
  fdr_threshold = fdr_threshold,
  number_association = FALSE,
  gene_set_name = gene_set_name,
  enrichr_set_name = enrichr_set_name,
  y_limits = c(0, 100)  # opcjonalnie – proporcje
)

p6 <- plot_AggregateEnrichr_summarySignatures_normalized(
  data = gene_summaryList_lung$n10_short_time,
  subdata = "summary_down",
  # plot_title = "Σ GR-related aggregated enrichments (Blood, normalized, down)",
  x_label = "",
  y_label = "Percent of genes from signature",
  palette_sources = c(
    "ChEA_2022" = "#08306b",
    "LINCS_L1000" = "#2171b5",
    "CellMarker_2024" = "#bdd7e7"
  ),
  axis_text_size = 20,
  axis_title_size = 20,
  legend_text_size = 18,
  signature_text_size = 20,
  fdr_threshold = fdr_threshold,
  number_association = FALSE,
  gene_set_name = gene_set_name,
  enrichr_set_name = enrichr_set_name,
  y_limits = c(0, 100)  # odpowiadające proporcjonalnemu zakresowi
)


(p1 + p2) / (p3 + p4) / (p5 + p6)


(p1 + p2) / (p3 + p4) / (p5 + p6) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")



# --- Dodanie pełnych opisów do lewej kolumny
p1 <- p1 + ggtitle("Signatures for neural tissue")
p3 <- p3 + ggtitle("Signatures for blood cells")
p5 <- p5 + ggtitle("Signatures for lung tissue")

# --- Połączenie wszystkich wykresów w układzie 3 wierszy, 2 kolumny
final_plot <- (p1 + p2) / (p3 + p4) / (p5 + p6) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    plot.title = element_text(
      size = 22,
      face = "plain",
      hjust = 0,        # wyrównanie do lewej
      margin = margin(b = 6)
    )
  )

final_plot
