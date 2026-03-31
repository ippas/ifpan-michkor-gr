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
    Organ = factor(Organ, levels = selected_organs)
  )

# opcjonalna kontrola duplikatów
duplicate_check <- df_plot %>%
  count(tissue, gene_ensembl, Organ) %>%
  filter(n > 1)

# =========================================================
# ANOVA
# =========================================================

anova_results <- df_plot %>%
  group_by(tissue) %>%
  anova_test(
    dv = log_expr,
    wid = gene_ensembl,
    within = Organ
  )

# =========================================================
# POSTHOC
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
    hf_dfd = as.numeric(hf_dfd)
  )

# =========================================================
# LABELS FOR GLOBAL ANOVA ON PLOT
# =========================================================

anova_plot_labels <- anova_results_summary_df %>%
  select(tissue, anova_p, anova_signif) %>%
  mutate(
    tissue = factor(tissue, levels = tissue_levels)
  )

label_y_df <- df_plot %>%
  group_by(tissue) %>%
  summarise(
    y_anova = max(log_expr, na.rm = TRUE) + 0.35,
    .groups = "drop"
  )

anova_plot_labels <- anova_plot_labels %>%
  left_join(label_y_df, by = "tissue") %>%
  mutate(
    x = 2
  )

# =========================================================
# POSTHOC LABELS FOR PLOT
# only neural, blood, lung
# only significant comparisons
# =========================================================

posthoc_plot_df <- posthoc_results %>%
  filter(tissue %in% c("neural", "blood", "lung")) %>%
  filter(p.adj < 0.05) %>%
  mutate(
    tissue = factor(tissue, levels = tissue_levels),
    p_label = p_to_stars(p.adj)
  ) %>%
  group_by(tissue) %>%
  arrange(group1, group2, .by_group = TRUE) %>%
  mutate(
    y.position = max(df_plot$log_expr[df_plot$tissue == unique(tissue)], na.rm = TRUE) +
      0.75 + seq(0, by = 0.40, length.out = n())
  ) %>%
  ungroup()

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
    size = 2
  ) +
  geom_text(
    data = anova_plot_labels,
    aes(x = x, y = y_anova, label = anova_signif),
    inherit.aes = FALSE,
    size = 6
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
  scale_color_manual(values = signature_colors) +
  scale_x_discrete(labels = organ_labels) +
  facet_wrap(
    ~ tissue,
    nrow = 1
  ) +
  coord_cartesian(
    ylim = c(
      min(df_plot$log_expr, na.rm = TRUE),
      max(
        c(
          anova_plot_labels$y_anova,
          if (nrow(posthoc_plot_df) > 0) posthoc_plot_df$y.position else NA_real_
        ),
        na.rm = TRUE
      ) + 0.25
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
posthoc_plot_df
plot_box