library(dplyr)
library(stringr)
library(purrr)
library(ggplot2)
library(rstatix)

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
    Organ = factor(
      Organ,
      levels = c(
        "brain",
        "bone marrow & lymphoid tissues",
        "respiratory system"
      )
    )
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

nova_results_df <-  tibble(
  tissue = anova_results$tissue,
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
    mauchly_df = map(
      mauchly_df,
      ~ .x %>%
        transmute(
          mauchly_effect = Effect,
          mauchly_W = W,
          mauchly_p = p,
          mauchly_signif = `p<.05`
        )
    ),
    sphericity_df = map(
      sphericity_df,
      ~ .x %>%
        transmute(
          correction_effect = Effect,
          gg_epsilon = GGe,
          gg_df = `DF[GG]`,
          gg_p = `p[GG]`,
          gg_signif = `p[GG]<.05`,
          hf_epsilon = HFe,
          hf_df = `DF[HF]`,
          hf_p = `p[HF]`,
          hf_signif = `p[HF]<.05`
        )
    )
  ) %>%
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
# PLOT
# =========================================================

plot_box <- ggplot(
  df_plot,
  aes(x = Organ, y = log_expr)
) +
  geom_boxplot(
    outlier.shape = NA,
    color = "black"
  ) +
  geom_jitter(
    aes(color = tissue),
    position = position_jitter(width = 0.18, height = 0),
    alpha = 0.8,
    size = 2
  ) +
  scale_color_manual(values = signature_colors) +
  scale_x_discrete(
    labels = c(
      "brain" = "Brain",
      "bone marrow & lymphoid tissues" = "Bone marrow\n& lymphoid tissues",
      "respiratory system" = "Respiratory system"
    )
  ) +
  facet_wrap(
    ~ tissue,
    nrow = 1
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
posthoc_results
plot_box