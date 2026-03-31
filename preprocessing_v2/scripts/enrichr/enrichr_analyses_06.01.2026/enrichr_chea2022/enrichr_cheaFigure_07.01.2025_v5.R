## ============================================================
## ADD: 3 dots per tile (brain / blood / lung) from GTEx medians
## - dot is filled if expressed (median_max > 1), hollow/grey if not
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})


top10_per_signature %>% 
  mutate(TF2 = TF) %>% 
  mutate(
    TF = case_when(
      TF == "LXR" ~ "NR1H3",
      TF == "SA1" ~ "STAG1",
      TF == "SMC1" ~ "SMC1A",
      TF == "CJUN" ~ "JUN",
      TF == "TCFCP2L1" ~ "TFCP2L1",
      TF == "Nerf2" ~ "NFE2L2",
      TF == "RING1B" ~ "RING1",
      TF == "AF4" ~ "AFF4",
      TRUE        ~ TF
    )
  ) -> top10_per_signature

top10_per_signature$TF

TF_expression_GTEx <- multi_downloadMedianSubtissuesGTEx_tissueSummaryExpression(
  gene_symbol = top10_per_signature$TF,
  verbose = TRUE,
  show_progress = TRUE,
  keep_unmapped = FALSE
)


## ----------------------------
## A) Prepare GTEx expression flags (TRUE/FALSE) per TF
## ----------------------------
TF_expr_flags <- TF_expression_GTEx %>%
  filter(tissue %in% c("Lung", "Brain", "Whole_Blood")) %>%
  mutate(
    tissue = recode(
      tissue,
      "Brain"       = "brain",
      "Lung"        = "lung",
      "Whole_Blood" = "blood"
    ),
    expressed = median_max > 1
  ) %>%
  group_by(query_gene_symbol, tissue) %>%
  summarise(expressed = any(expressed), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from  = tissue,
    values_from = expressed,
    values_fill = FALSE
  )

## ----------------------------
## B) Join to heatmap DF and create dot coordinates (3 per tile)
## ----------------------------
dot_pal <- c(
  brain = "#2E86FF",
  lung  = "#00A676",
  blood = "#E31A1C"
)

dots_df <- hm_df %>%
  left_join(TF_expr_flags, by = c("TF" = "query_gene_symbol")) %>%
  mutate(
    brain = ifelse(is.na(brain), FALSE, brain),
    lung  = ifelse(is.na(lung),  FALSE, lung),
    blood = ifelse(is.na(blood), FALSE, blood)
  ) %>%
  tidyr::pivot_longer(
    cols = c(brain, lung, blood),
    names_to = "tissue",
    values_to = "expressed"
  ) %>%
  mutate(
    tissue = factor(tissue, levels = c("brain", "lung", "blood")),
    ## numeric position of tile center:
    x_center = x_num,
    y_center = y_num,
    ## place 3 dots as a small row near the TOP of tile
    ## (remember scale_y_reverse(): smaller y is higher on screen)
    x_dot = x_center + c(brain = -0.23, lung = 0.00, blood = 0.23)[as.character(tissue)],
    y_dot = y_center - 0.30
  )

## ----------------------------
## C) Add layers to your existing plot object
## ----------------------------
p_heat_cat_dots <- p_heat_cat +
  ## hollow dot as "not expressed" background (so it’s visible on dark bins too)
  geom_point(
    data = dots_df,
    inherit.aes = FALSE,
    aes(x = x_dot, y = y_dot),
    shape = 21,
    size  = 2.2,
    stroke = 0.7,
    fill  = "white",
    color = "grey35",
    na.rm = TRUE
  ) +
  ## filled colored dot where expressed == TRUE
  geom_point(
    data = dots_df %>% filter(expressed),
    inherit.aes = FALSE,
    aes(x = x_dot, y = y_dot, fill = tissue),
    shape  = 21,
    size   = 2.2,
    stroke = 0.7,
    color  = "grey10",
    na.rm  = TRUE
  ) +
  scale_fill_manual(values = c(fill_pal, dot_pal), breaks = names(dot_pal), guide = "none")
# ^ ważne: nie chcemy mieszać legendy FDR z dotami, więc guide = "none"

p_heat_cat_dots
