library(dplyr)
library(stringr)
library(purrr)
library(tidyr)
library(ggplot2)
library(rstatix)
library(tibble)
library(ggpubr)

# =========================================================
# SETTINGS
# =========================================================

signature_colors <- c(
  systemic = "#335C67",
  neural   = "#E09F3E",
  blood    = "#9E2A2B",
  lung     = "#540B0E"
)

selected_organs <- c(
  "brain",
  "bone marrow & lymphoid tissues",
  "respiratory system"
)

tissue_levels <- c(
  "systemic",
  "neural",
  "blood",
  "lung"
)

organ_labels <- c(
  "brain" = "Brain",
  "bone marrow & lymphoid tissues" = "Bone marrow\n& lymphoid tissues",
  "respiratory system" = "Respiratory system"
)

# =========================================================
# HELPERS
# =========================================================

get_overlap_genes <- function(signature_name) {
  cleanOverlapResults_grSystemicTissues_ICD10F %>%
    filter(grepl(signature_name, grSignature)) %>%
    filter(p_value < 0.05) %>%
    filter(observed_overlap >= 3) %>%
    pull(overlap_genes) %>%
    strsplit(",") %>%
    unlist() %>%
    trimws() %>%
    unique()
}

prepare_hpa_df <- function(gene_list, tissue_name) {
  HumanProteinAtlas_rnaTissueConsensus %>%
    filter(Organ %in% selected_organs) %>%
    select(Gene, gene_symbol, Organ, organ_nTPM_max) %>%
    filter(gene_symbol %in% gene_list) %>%
    distinct() %>%
    rename(gene_ensembl = Gene) %>%
    mutate(
      log_expr = log2(organ_nTPM_max + 1),
      tissue   = tissue_name
    )
}

p_to_stars <- function(p) {
  case_when(
    is.na(p)  ~ "n.s.",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "n.s."
  )
}

# =========================================================
# GENE LISTS
# =========================================================

genes_systemic <- get_overlap_genes("systemic")
genes_neural   <- get_overlap_genes("neural")
genes_blood    <- get_overlap_genes("blood")
genes_lung     <- get_overlap_genes("lung")

# =========================================================
# DATA FOR PLOT AND STATISTICS
# =========================================================

df_plot <- bind_rows(
  prepare_hpa_df(genes_systemic, "systemic"),
  prepare_hpa_df(genes_neural,   "neural"),
  prepare_hpa_df(genes_blood,    "blood"),
  prepare_hpa_df(genes_lung,     "lung")
) %>%
  distinct(tissue, gene_ensembl, gene_symbol, Organ, .keep_all = TRUE) %>%
  group_by(tissue, gene_ensembl) %>%
  filter(n_distinct(Organ) == 3) %>%
  ungroup() %>%
  mutate(
    tissue = factor(tissue, levels = tissue_levels),
    Organ  = factor(Organ, levels = selected_organs)
  )

# opcjonalna kontrola
duplicate_check <- df_plot %>%
  count(tissue, gene_ensembl, Organ) %>%
  filter(n > 1)

# =========================================================
# REPEATED-MEASURES ANOVA
# =========================================================

anova_results <- df_plot %>%
  group_by(tissue) %>%
  anova_test(
    dv = log_expr,
    wid = gene_ensembl,
    within = Organ
  )

# =========================================================
# POST HOC
# =========================================================

posthoc_results <- df_plot %>%
  group_by(tissue) %>%
  pairwise_t_test(
    log_expr ~ Organ,
    paired = TRUE,
    p.adjust.method = "BH"
  )

# =========================================================
# FORMATTED ANOVA SUMMARY TABLE
# =========================================================

anova_results_summary_df <- tibble(
  tissue = anova_results$tissue,
  anova_df = map(
    anova_results$anova,
    ~ as_tibble(as.data.frame(.x[["ANOVA"]]))
  ),
  mauchly_df = map(
    anova_results$anova,
    ~ as_tibble(as.data.frame(.x[["Mauchly's Test for Sphericity"]]))
  ),
  sphericity_df = map(
    anova_results$anova,
    ~ as_tibble(as.data.frame(.x[["Sphericity Corrections"]]))
  )
) %>%
  mutate(
    anova_df = map(
      anova_df,
      ~ .x %>%
        transmute(
          anova_effect = Effect,
          anova_dfn    = DFn,
          anova_dfd    = DFd,
          anova_F      = F,
          anova_p      = p,
          anova_signif = p_to_stars(p),
          ges          = ges
        )
    ),
    mauchly_df = map(
      mauchly_df,
      ~ .x %>%
        transmute(
          mauchly_W = W,
          mauchly_p = p
        )
    ),
    sphericity_df = map(
      sphericity_df,
      ~ .x %>%
        transmute(
          gg_epsilon = GGe,
          gg_df      = `DF[GG]`,
          gg_p       = `p[GG]`,
          gg_signif  = p_to_stars(`p[GG]`),
          hf_epsilon = HFe,
          hf_df      = `DF[HF]`,
          hf_p       = `p[HF]`,
          hf_signif  = p_to_stars(`p[HF]`)
        )
    )
  ) %>%
  unnest(anova_df) %>%
  unnest(mauchly_df) %>%
  unnest(sphericity_df) %>%
  separate(gg_df, into = c("gg_dfn", "gg_dfd"), sep = ",\\s*") %>%
  separate(hf_df, into = c("hf_dfn", "hf_dfd"), sep = ",\\s*") %>%
  mutate(
    gg_dfn = as.numeric(gg_dfn),
    gg_dfd = as.numeric(gg_dfd),
    hf_dfn = as.numeric(hf_dfn),
    hf_dfd = as.numeric(hf_dfd),
    tissue = factor(tissue, levels = tissue_levels)
  )

# =========================================================
# PANEL MAX VALUES
# =========================================================

panel_max_df <- df_plot %>%
  group_by(tissue) %>%
  summarise(
    panel_y_max = max(log_expr, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    tissue = factor(tissue, levels = tissue_levels)
  )

# =========================================================
# GLOBAL ANOVA LABELS
# =========================================================

anova_plot_labels <- anova_results_summary_df %>%
  select(tissue, anova_p, anova_signif) %>%
  distinct() %>%
  left_join(panel_max_df, by = "tissue") %>%
  mutate(
    x = 2
  )

# =========================================================
# POST HOC LABELS
# show only where global ANOVA is significant
# but include ALL pairwise comparisons, also n.s.
# =========================================================

significant_anova_tissues <- anova_plot_labels %>%
  filter(anova_p < 0.05) %>%
  pull(tissue) %>%
  as.character()

posthoc_plot_df <- posthoc_results %>%
  filter(tissue %in% significant_anova_tissues) %>%
  mutate(
    tissue = factor(tissue, levels = tissue_levels),
    p_label = p_to_stars(p.adj)
  ) %>%
  left_join(panel_max_df, by = "tissue") %>%
  group_by(tissue) %>%
  arrange(group1, group2, .by_group = TRUE) %>%
  mutate(
    y.position = panel_y_max + 0.35 + seq(0, by = 0.40, length.out = n())
  ) %>%
  ungroup()

# =========================================================
# MOVE ANOVA ABOVE POST HOC
# =========================================================

anova_y_df <- posthoc_plot_df %>%
  group_by(tissue) %>%
  summarise(
    y_anova = max(y.position) + 0.40,
    .groups = "drop"
  )

anova_plot_labels <- anova_plot_labels %>%
  left_join(anova_y_df, by = "tissue") %>%
  mutate(
    y_anova = ifelse(
      is.na(y_anova),
      panel_y_max + 0.35,
      y_anova
    )
  )

# =========================================================
# Y LIMIT
# =========================================================

y_upper_limit <- max(
  c(
    anova_plot_labels$y_anova,
    if (nrow(posthoc_plot_df) > 0) posthoc_plot_df$y.position else NA_real_
  ),
  na.rm = TRUE
) + 0.25

# =========================================================
# FINAL PLOT
# =========================================================

plot_box <- ggplot(
  df_plot,
  aes(x = Organ, y = log_expr)
) +
  geom_boxplot(
    outlier.shape = NA,
    color = "black",
    width = 0.65
  ) +
  geom_jitter(
    aes(color = tissue),
    position = position_jitter(width = 0.18, height = 0),
    alpha = 0.8,
    size = 1
  ) +
  stat_pvalue_manual(
    posthoc_plot_df,
    label = "p_label",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    bracket.size = 0.5,
    size = 5,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = anova_plot_labels,
    aes(x = x, y = y_anova, label = anova_signif),
    inherit.aes = FALSE,
    size = 6
  ) +
  scale_color_manual(values = signature_colors) +
  scale_x_discrete(labels = organ_labels) +
  facet_wrap(
    ~ tissue,
    nrow = 1
  ) +
  coord_cartesian(
    ylim = c(
      min(df_plot$log_expr, na.rm = TRUE),
      y_upper_limit
    )
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  ) +
  labs(
    x = "Organ",
    y = "log2(max nTPM + 1)",
    title = "Expression of overlap genes across organs"
  )

# =========================================================
# OUTPUT OBJECTS
# =========================================================

df_plot
anova_results
anova_results_summary_df
posthoc_results
anova_plot_labels
posthoc_plot_df
plot_box


# ##############################################################################
# SAVE TO SVG
# ##############################################################################
dev.off()

svg(
  filename = "/home/mateusz/projects/ifpan-michkor-gr/results_v2/HumanProteinAtlas/boxplotHPAexpression_overlapResults_grSystemicTissues_v1_12.03.2026/boxplotHPAexpression_overlapResults_grSystemicTissues_v1_12.03.2026.svg",
  width = 8.5,
  height = 5.5
)

# tutaj wstawiasz kod tworzący wykres
print(plot_box)

dev.off()


# ##############################################################################
# SAVE TO XLSX
# ##############################################################################
library(dplyr)
library(tidyr)
library(purrr)
library(openxlsx)

# =========================================================
# PATHS
# =========================================================

out_dir <- "/home/mateusz/projects/ifpan-michkor-gr/results_v2/HumanProteinAtlas/boxplotHPAexpression_overlapResults_grSystemicTissues_v1_12.03.2026"

out_file <- file.path(
  out_dir,
  "supplement_HPAexpressionAndStatistics_overlapResultsGrSystemicTissues_v1_12.03.2026.xlsx"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# PREPARE SUPPLEMENT TABLE FROM df_plot
# =========================================================

supplement_df <- df_plot %>%
  rename(
    ensembl_id = gene_ensembl,
    signatureType = tissue,
    nTPM_max = organ_nTPM_max
  ) %>%
  mutate(
    Organ_clean = case_when(
      Organ == "brain" ~ "brain",
      Organ == "bone marrow & lymphoid tissues" ~ "bone_marrow_lymphoid_tissues",
      Organ == "respiratory system" ~ "respiratory_system"
    )
  ) %>%
  select(
    ensembl_id,
    gene_symbol,
    signatureType,
    Organ_clean,
    nTPM_max,
    log_expr
  ) %>%
  pivot_wider(
    id_cols = c(ensembl_id, gene_symbol, signatureType),
    names_from = Organ_clean,
    values_from = c(nTPM_max, log_expr),
    names_glue = "{Organ_clean}_{.value}"
  ) %>%
  select(
    ensembl_id,
    gene_symbol,
    signatureType,
    brain_nTPM_max,
    bone_marrow_lymphoid_tissues_nTPM_max,
    respiratory_system_nTPM_max,
    brain_log_expr,
    bone_marrow_lymphoid_tissues_log_expr,
    respiratory_system_log_expr
  ) %>%
  mutate(
    signatureType = factor(
      as.character(signatureType),
      levels = c("systemic", "neural", "blood", "lung")
    )
  ) %>%
  arrange(signatureType, gene_symbol, ensembl_id)

supplement_df_list <- supplement_df %>%
  split(.$signatureType) %>%
  .[c("systemic", "neural", "blood", "lung")] %>%
  map(~ .x %>% select(-signatureType))

# =========================================================
# PREPARE STATISTICS TABLES
# =========================================================

anova_results_df <- anova_results_summary_df %>%
  mutate(
    tissue = factor(as.character(tissue), levels = c("systemic", "neural", "blood", "lung"))
  ) %>%
  arrange(tissue)

posthoc_results_df <- posthoc_results %>%
  mutate(
    tissue = factor(as.character(tissue), levels = c("systemic", "neural", "blood", "lung"))
  ) %>%
  arrange(tissue, group1, group2)

# =========================================================
# WORKBOOK AND STYLES
# =========================================================

wb <- createWorkbook()

header_style <- createStyle(
  textDecoration = "bold",
  halign = "center",
  valign = "center"
)

title_style <- createStyle(
  textDecoration = "bold",
  fontSize = 12
)

log_style <- createStyle(numFmt = "0.000")

# =========================================================
# 1. STATISTICS SHEET
# =========================================================

addWorksheet(wb, "STATISTICS")

writeData(
  wb,
  sheet = "STATISTICS",
  x = "anova_results_summary_df",
  startRow = 1,
  startCol = 1,
  colNames = FALSE
)

addStyle(
  wb,
  sheet = "STATISTICS",
  style = title_style,
  rows = 1,
  cols = 1,
  gridExpand = FALSE
)

writeData(
  wb,
  sheet = "STATISTICS",
  x = anova_results_df,
  startRow = 2,
  startCol = 1,
  withFilter = FALSE
)

addStyle(
  wb,
  sheet = "STATISTICS",
  style = header_style,
  rows = 2,
  cols = 1:ncol(anova_results_df),
  gridExpand = TRUE
)

writeData(
  wb,
  sheet = "STATISTICS",
  x = "posthoc_results",
  startRow = 10,
  startCol = 1,
  colNames = FALSE
)

addStyle(
  wb,
  sheet = "STATISTICS",
  style = title_style,
  rows = 10,
  cols = 1,
  gridExpand = FALSE
)

writeData(
  wb,
  sheet = "STATISTICS",
  x = posthoc_results_df,
  startRow = 11,
  startCol = 1,
  withFilter = FALSE
)

addStyle(
  wb,
  sheet = "STATISTICS",
  style = header_style,
  rows = 11,
  cols = 1:ncol(posthoc_results_df),
  gridExpand = TRUE
)

freezePane(
  wb,
  sheet = "STATISTICS",
  firstRow = TRUE
)

setColWidths(
  wb,
  sheet = "STATISTICS",
  cols = 1:max(ncol(anova_results_df), ncol(posthoc_results_df)),
  widths = "auto"
)

# =========================================================
# 2. SIGNATURE SHEETS
# =========================================================

for (sheet_name in c("systemic", "neural", "blood", "lung")) {
  
  df <- supplement_df_list[[sheet_name]]
  
  addWorksheet(wb, sheet_name)
  
  writeData(
    wb,
    sheet = sheet_name,
    x = df,
    startRow = 1,
    startCol = 1,
    withFilter = FALSE
  )
  
  addStyle(
    wb,
    sheet = sheet_name,
    style = header_style,
    rows = 1,
    cols = 1:ncol(df),
    gridExpand = TRUE
  )
  
  addFilter(
    wb,
    sheet = sheet_name,
    rows = 1,
    cols = 1:ncol(df)
  )
  
  freezePane(
    wb,
    sheet = sheet_name,
    firstRow = TRUE
  )
  
  setColWidths(
    wb,
    sheet = sheet_name,
    cols = 1:ncol(df),
    widths = "auto"
  )
  
  log_cols <- grep("log_expr$", names(df))
  
  addStyle(
    wb,
    sheet = sheet_name,
    style = log_style,
    rows = 2:(nrow(df) + 1),
    cols = log_cols,
    gridExpand = TRUE
  )
}

# =========================================================
# SAVE
# =========================================================

saveWorkbook(wb, out_file, overwrite = TRUE)

out_file
