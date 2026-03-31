# ##############################################################################
# ---- min 3 publication ----
# ##############################################################################

# ##############################################################################
# ---- blood ----
# ##############################################################################
p1 <- plot_AggregateEnrichr_summarySignatures(
  data = gene_summaryList_blood$n10_short_time,
  subdata = "summary_up",
  # plot_title = "Σ GR-related aggregated enrichments (Blood, n10_short_time)",
  x_label = "",
  y_label = "Number of unique genes",
  palette_sources = c(
    "ChEA_2022" = "#4b0b0b",
    "LINCS_L1000" = "#9c5757",
    "CellMarker_2024" = "#f0c1c1"
  ),
  axis_text_size = 20,
  axis_title_size = 20,
  legend_text_size = 18,
  signature_text_size = 20,
  number_association = FALSE
)


p2 <- plot_AggregateEnrichr_summarySignatures(
  data = gene_summaryList_blood$n10_short_time,
  subdata = "summary_down",
  # plot_title = "Σ GR-related aggregated enrichments (Blood, n10_short_time)",
  x_label = "",
  y_label = "Number of unique genes",
  palette_sources = c(
    "ChEA_2022" = "#08306b",
    "LINCS_L1000" = "#2171b5",
    "CellMarker_2024" = "#bdd7e7"
  ),
  axis_text_size = 20,
  axis_title_size = 20,
  legend_text_size = 18,
  signature_text_size = 20,
  y_limits = c(0, 40),
  number_association = FALSE
)

p1 + p2
# ##############################################################################
# ---- brain ----
# ##############################################################################
p1 <- plot_AggregateEnrichr_summarySignatures(
  data = gene_summaryList_brain$n10_short_time,
  subdata = "summary_up",
  # plot_title = "Σ GR-related aggregated enrichments (Blood, n10_short_time)",
  x_label = "",
  y_label = "Number of unique genes",
  palette_sources = c(
    "ChEA_2022" = "#4b0b0b",
    "LINCS_L1000" = "#9c5757",
    "CellMarker_2024" = "#f0c1c1"
  ),
  axis_text_size = 20,
  axis_title_size = 20,
  legend_text_size = 18,
  signature_text_size = 20,
  number_association = FALSE
)


p2 <- plot_AggregateEnrichr_summarySignatures(
  data = gene_summaryList_brain$n10_short_time,
  subdata = "summary_down",
  # plot_title = "Σ GR-related aggregated enrichments (Blood, n10_short_time)",
  x_label = "",
  y_label = "Number of unique genes",
  palette_sources = c(
    "ChEA_2022" = "#08306b",
    "LINCS_L1000" = "#2171b5",
    "CellMarker_2024" = "#bdd7e7"
  ),
  axis_text_size = 20,
  axis_title_size = 20,
  legend_text_size = 18,
  signature_text_size = 20,
  number_association = FALSE
)


p1 + p2
# ##############################################################################
# ---- lung ----
# ##############################################################################
p1 <- plot_AggregateEnrichr_summarySignatures(
  data = gene_summaryList_lung$n10_short_time,
  subdata = "summary_up",
  # plot_title = "Σ GR-related aggregated enrichments (Blood, n10_short_time)",
  x_label = "",
  y_label = "Number of unique genes",
  palette_sources = c(
    "ChEA_2022" = "#4b0b0b",
    "LINCS_L1000" = "#9c5757",
    "CellMarker_2024" = "#f0c1c1"
  ),
  axis_text_size = 20,
  axis_title_size = 20,
  legend_text_size = 18,
  signature_text_size = 20,
  number_association = FALSE
)


p2 <- plot_AggregateEnrichr_summarySignatures(
  data = gene_summaryList_lung$n10_short_time,
  subdata = "summary_down",
  # plot_title = "Σ GR-related aggregated enrichments (Blood, n10_short_time)",
  x_label = "",
  y_label = "Number of unique genes",
  palette_sources = c(
    "ChEA_2022" = "#08306b",
    "LINCS_L1000" = "#2171b5",
    "CellMarker_2024" = "#bdd7e7"
  ),
  axis_text_size = 20,
  axis_title_size = 20,
  legend_text_size = 18,
  signature_text_size = 20,
  number_association = FALSE
)

p1 + p2


# ##############################################################################
# ---- min 4 publication ----
# ##############################################################################

# ##############################################################################
# ---- blood ----
# ##############################################################################
p1 <- plot_AggregateEnrichr_summarySignatures(
  data = gene_summaryList_blood$n10_short_time,
  subdata = "summary_up",
  # plot_title = "Σ GR-related aggregated enrichments (Blood, n10_short_time)",
  x_label = "",
  y_label = "Number of unique genes",
  palette_sources = c(
    "ChEA_2022" = "#4b0b0b",
    "LINCS_L1000" = "#9c5757",
    "CellMarker_2024" = "#f0c1c1"
  ),
  axis_text_size = 20,
  axis_title_size = 20,
  legend_text_size = 18,
  signature_text_size = 20,
  number_association = FALSE,
  y_limits = c(0, 25),
  gene_set_name = "genes_4pub",
  enrichr_set_name = "enrichr_4pub"
)


p2 <- plot_AggregateEnrichr_summarySignatures(
  data = gene_summaryList_blood$n10_short_time,
  subdata = "summary_down",
  # plot_title = "Σ GR-related aggregated enrichments (Blood, n10_short_time)",
  x_label = "",
  y_label = "Number of unique genes",
  palette_sources = c(
    "ChEA_2022" = "#08306b",
    "LINCS_L1000" = "#2171b5",
    "CellMarker_2024" = "#bdd7e7"
  ),
  axis_text_size = 20,
  axis_title_size = 20,
  legend_text_size = 18,
  signature_text_size = 20,
  y_limits = c(0, 25),
  number_association = FALSE,
  gene_set_name = "genes_4pub",
  enrichr_set_name = "enrichr_4pub"
)

p1 + p2
# ##############################################################################
# ---- brain ----
# ##############################################################################
p1 <- plot_AggregateEnrichr_summarySignatures(
  data = gene_summaryList_brain$n10_short_time,
  subdata = "summary_up",
  # plot_title = "Σ GR-related aggregated enrichments (Blood, n10_short_time)",
  x_label = "",
  y_label = "Number of unique genes",
  palette_sources = c(
    "ChEA_2022" = "#4b0b0b",
    "LINCS_L1000" = "#9c5757",
    "CellMarker_2024" = "#f0c1c1"
  ),
  axis_text_size = 20,
  axis_title_size = 20,
  legend_text_size = 18,
  signature_text_size = 20,
  y_limits = c(0, 50),
  number_association = FALSE,
  gene_set_name = "genes_4pub",
  enrichr_set_name = "enrichr_4pub"
)


p2 <- plot_AggregateEnrichr_summarySignatures(
  data = gene_summaryList_brain$n10_short_time,
  subdata = "summary_down",
  # plot_title = "Σ GR-related aggregated enrichments (Blood, n10_short_time)",
  x_label = "",
  y_label = "Number of unique genes",
  palette_sources = c(
    "ChEA_2022" = "#08306b",
    "LINCS_L1000" = "#2171b5",
    "CellMarker_2024" = "#bdd7e7"
  ),
  axis_text_size = 20,
  axis_title_size = 20,
  legend_text_size = 18,
  signature_text_size = 20,
  y_limits = c(0, 50),
  number_association = FALSE,
  gene_set_name = "genes_4pub",
  enrichr_set_name = "enrichr_4pub"
)


p1 + p2
# ##############################################################################
# ---- lung ----
# ##############################################################################
p1 <- plot_AggregateEnrichr_summarySignatures(
  data = gene_summaryList_lung$n10_short_time,
  subdata = "summary_up",
  # plot_title = "Σ GR-related aggregated enrichments (Blood, n10_short_time)",
  x_label = "",
  y_label = "Number of unique genes",
  palette_sources = c(
    "ChEA_2022" = "#4b0b0b",
    "LINCS_L1000" = "#9c5757",
    "CellMarker_2024" = "#f0c1c1"
  ),
  axis_text_size = 20,
  axis_title_size = 20,
  legend_text_size = 18,
  signature_text_size = 20,
  y_limits = c(0, 100),
  number_association = FALSE,
  gene_set_name = "genes_4pub",
  enrichr_set_name = "enrichr_4pub"
)


p2 <- plot_AggregateEnrichr_summarySignatures(
  data = gene_summaryList_lung$n10_short_time,
  subdata = "summary_down",
  # plot_title = "Σ GR-related aggregated enrichments (Blood, n10_short_time)",
  x_label = "",
  y_label = "Number of unique genes",
  palette_sources = c(
    "ChEA_2022" = "#08306b",
    "LINCS_L1000" = "#2171b5",
    "CellMarker_2024" = "#bdd7e7"
  ),
  axis_text_size = 20,
  axis_title_size = 20,
  legend_text_size = 18,
  signature_text_size = 20,
  y_limits = c(0, 100),
  number_association = FALSE,
  gene_set_name = "genes_4pub",
  enrichr_set_name = "enrichr_4pub"
)

p1 + p2


# ##############################################################################
# ---- min 5 publication ----
# ##############################################################################

# ##############################################################################
# ---- blood ----
# ##############################################################################
p1 <- plot_AggregateEnrichr_summarySignatures(
  data = gene_summaryList_blood$n10_short_time,
  subdata = "summary_up",
  # plot_title = "Σ GR-related aggregated enrichments (Blood, n10_short_time)",
  x_label = "",
  y_label = "Number of unique genes",
  palette_sources = c(
    "ChEA_2022" = "#4b0b0b",
    "LINCS_L1000" = "#9c5757",
    "CellMarker_2024" = "#f0c1c1"
  ),
  axis_text_size = 20,
  axis_title_size = 20,
  legend_text_size = 18,
  signature_text_size = 20,
  number_association = FALSE,
  y_limits = c(0, 6),
  gene_set_name = "genes_5pub",
  enrichr_set_name = "enrichr_5pub"
)


p2 <- plot_AggregateEnrichr_summarySignatures(
  data = gene_summaryList_blood$n10_short_time,
  subdata = "summary_down",
  # plot_title = "Σ GR-related aggregated enrichments (Blood, n10_short_time)",
  x_label = "",
  y_label = "Number of unique genes",
  palette_sources = c(
    "ChEA_2022" = "#08306b",
    "LINCS_L1000" = "#2171b5",
    "CellMarker_2024" = "#bdd7e7"
  ),
  axis_text_size = 20,
  axis_title_size = 20,
  legend_text_size = 18,
  signature_text_size = 20,
  y_limits = c(0, 6),
  number_association = FALSE,
  gene_set_name = "genes_5pub",
  enrichr_set_name = "enrichr_5pub"
)

p1 + p2
# ##############################################################################
# ---- brain ----
# ##############################################################################
p1 <- plot_AggregateEnrichr_summarySignatures(
  data = gene_summaryList_brain$n10_short_time,
  subdata = "summary_up",
  # plot_title = "Σ GR-related aggregated enrichments (Blood, n10_short_time)",
  x_label = "",
  y_label = "Number of unique genes",
  palette_sources = c(
    "ChEA_2022" = "#4b0b0b",
    "LINCS_L1000" = "#9c5757",
    "CellMarker_2024" = "#f0c1c1"
  ),
  axis_text_size = 20,
  axis_title_size = 20,
  legend_text_size = 18,
  signature_text_size = 20,
  y_limits = c(0, 10),
  number_association = FALSE,
  gene_set_name = "genes_5pub",
  enrichr_set_name = "enrichr_5pub"
)


p2 <- plot_AggregateEnrichr_summarySignatures(
  data = gene_summaryList_brain$n10_short_time,
  subdata = "summary_down",
  # plot_title = "Σ GR-related aggregated enrichments (Blood, n10_short_time)",
  x_label = "",
  y_label = "Number of unique genes",
  palette_sources = c(
    "ChEA_2022" = "#08306b",
    "LINCS_L1000" = "#2171b5",
    "CellMarker_2024" = "#bdd7e7"
  ),
  axis_text_size = 20,
  axis_title_size = 20,
  legend_text_size = 18,
  signature_text_size = 20,
  y_limits = c(0, 6),
  number_association = FALSE,
  gene_set_name = "genes_5pub",
  enrichr_set_name = "enrichr_5pub"
)


p1 + p2
# ##############################################################################
# ---- lung ----
# ##############################################################################
p1 <- plot_AggregateEnrichr_summarySignatures(
  data = gene_summaryList_lung$n10_short_time,
  subdata = "summary_up",
  # plot_title = "Σ GR-related aggregated enrichments (Blood, n10_short_time)",
  x_label = "",
  y_label = "Number of unique genes",
  palette_sources = c(
    "ChEA_2022" = "#4b0b0b",
    "LINCS_L1000" = "#9c5757",
    "CellMarker_2024" = "#f0c1c1"
  ),
  axis_text_size = 20,
  axis_title_size = 20,
  legend_text_size = 18,
  signature_text_size = 20,
  y_limits = c(0, 35),
  number_association = FALSE,
  gene_set_name = "genes_5pub",
  enrichr_set_name = "enrichr_5pub"
)


p2 <- plot_AggregateEnrichr_summarySignatures(
  data = gene_summaryList_lung$n10_short_time,
  subdata = "summary_down",
  # plot_title = "Σ GR-related aggregated enrichments (Blood, n10_short_time)",
  x_label = "",
  y_label = "Number of unique genes",
  palette_sources = c(
    "ChEA_2022" = "#08306b",
    "LINCS_L1000" = "#2171b5",
    "CellMarker_2024" = "#bdd7e7"
  ),
  axis_text_size = 20,
  axis_title_size = 20,
  legend_text_size = 18,
  signature_text_size = 20,
  y_limits = c(0, 35),
  number_association = FALSE,
  gene_set_name = "genes_5pub",
  enrichr_set_name = "enrichr_5pub"
)

p1 + p2

