## ============================================================
## TF heatmap (DISCRETE FDR BINS) + nuclear receptor bottom strip
## + 3 GTEx expression dots per tile (Brain / Lung / Blood)
## - bins include 1e-6
## - nuclear receptor list: nuclearReceptors_GO0004879
## - purple strip is at the BOTTOM of each tile (with y reversed)
## - dots: hollow if NOT expressed; colored if expressed (median_max > 1)
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
})


## ----------------------------
## 3) Build one table: top10 per signature
## ----------------------------
top10_per_signature <- res$enrichr$overlap_only %>%
  imap_dfr(~ .x %>%
             mutate(
               TF  = str_extract(Term, "^[^ ]+"),
               FDR = as.numeric(FDR)
             ) %>%
             group_by(TF) %>%
             slice_min(order_by = FDR, n = 1, with_ties = FALSE) %>%
             ungroup() %>%
             arrange(FDR) %>%
             slice_head(n = 10) %>%
             mutate(rank = row_number()) %>%
             mutate(signature = .y, .before = rank) %>%
             select(signature, rank, TF, everything())
  )


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
      TF == "NRF2" ~ "NFE2L2",
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
## 0) REQUIRED INPUTS
## ----------------------------
## Make sure these exist in your environment:
## - res$enrichr$overlap_only  (list: signature -> enrichr table)
## - nuclearReceptors_GO0004879 (character vector of TFs)
## - TF_expression_GTEx (dataframe with columns: TF, tissue, median_max)
##   NOTE: if TF column has a different name in your GTEx df, rename it to TF.

## ----------------------------
## 1) Signature order
## ----------------------------
signature_order <- c(
  "global_up", "global_down",
  "brain_up",  "brain_down",
  "blood_up",  "blood_down",
  "lung_up",   "lung_down"
)

## ----------------------------
## 2) FDR bins (include 1e-6)
## ----------------------------
fdr_breaks <- c(0, 1e-30, 1e-20, 1e-10, 1e-6, 1e-2, Inf)

fdr_labels <- c(
  "≤1e-30",
  "(1e-30, 1e-20]",
  "(1e-20, 1e-10]",
  "(1e-10, 1e-6]",
  "(1e-6, 0.01]",
  ">0.01"
)

fill_pal <- c(
  ">0.01"          = "white",
  "(1e-6, 0.01]"   = "#fae8e6",
  "(1e-10, 1e-6]"  = "#eba8a1",
  "(1e-20, 1e-10]" = "#d46b66",
  "(1e-30, 1e-20]" = "#b01728",
  "≤1e-30"         = "#5a0011"
)

text_pal <- c(
  ">0.01"          = "black",
  "(1e-6, 0.01]"   = "black",
  "(1e-10, 1e-6]"  = "black",
  "(1e-20, 1e-10]" = "black",
  "(1e-30, 1e-20]" = "white",
  "≤1e-30"         = "white"
)


## ----------------------------
## 4) Heatmap dataframe (complete grid + bins + coords)
## ----------------------------
hm_df <- top10_per_signature %>%
  mutate(
    signature = factor(signature, levels = signature_order),
    rank      = factor(rank, levels = 1:10, ordered = TRUE)
  ) %>%
  complete(signature, rank) %>%
  mutate(
    signature = factor(signature, levels = signature_order),
    rank      = factor(rank, levels = 1:10, ordered = TRUE),
    FDR_num   = as.numeric(FDR),
    fdr_bin   = cut(
      FDR_num,
      breaks = fdr_breaks,
      labels = fdr_labels,
      include.lowest = TRUE,
      right = TRUE
    ),
    fdr_bin   = factor(fdr_bin, levels = rev(fdr_labels)),
    text_col  = text_pal[as.character(fdr_bin)],
    is_NR     = !is.na(TF) & (TF %in% nuclearReceptors_GO0004879),
    x_num     = as.integer(signature),
    y_num     = as.integer(rank)
  )

## ----------------------------
## 5) GTEx expression flags per TF (Brain/Lung/Blood)
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
  pivot_wider(
    names_from  = tissue,
    values_from = expressed,
    values_fill = FALSE
  )

## ----------------------------
## 6) Build dots df: 3 dots per tile
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
  pivot_longer(
    cols = c(brain, lung, blood),
    names_to = "tissue",
    values_to = "expressed"
  ) %>%
  mutate(
    tissue = factor(tissue, levels = c("brain", "lung", "blood")),
    x_center = x_num,
    y_center = y_num,
    ## dot positions (small row near top of tile)
    x_dot = x_center + c(brain = -0.23, lung = 0.00, blood = 0.23)[as.character(tissue)],
    y_dot = y_center - 0.30
  )

## ----------------------------
## 7) Plot
## ----------------------------
nr_purple <- "#B388FF"  # lighter violet
strip_h   <- 0.10

dot_size  <- 3.0
dot_stroke <- 0.8

p_heat_cat <- ggplot(hm_df, aes(x = signature, y = y_num, fill = fdr_bin)) +
  ## base tiles
  geom_tile(color = "white", linewidth = 0.6, na.rm = FALSE) +
  
  ## nuclear receptor strip at the BOTTOM of the visible tile
  ## NOTE: because we use scale_y_reverse(), "bottom on screen" = higher y-values.
  geom_rect(
    data = hm_df %>% filter(is_NR),
    inherit.aes = FALSE,
    aes(
      xmin = x_num - 0.5,
      xmax = x_num + 0.5,
      ymin = y_num + 0.5 - strip_h,
      ymax = y_num + 0.5
    ),
    fill = nr_purple,
    color = NA
  ) +
  
  ## GTEx dots (draw hollow first, then filled for expressed)
  geom_point(
    data = dots_df,
    inherit.aes = FALSE,
    aes(x = x_dot, y = y_dot),
    shape  = 21,
    size   = dot_size,
    stroke = dot_stroke,
    fill   = "white",
    color  = "grey40",
    na.rm  = TRUE
  ) +
  geom_point(
    data = dots_df %>% filter(expressed),
    inherit.aes = FALSE,
    aes(x = x_dot, y = y_dot, color = tissue),
    shape  = 16,
    size   = dot_size,
    na.rm  = TRUE
  ) +
  
  ## TF labels
  geom_text(
    aes(label = TF, color = text_col),
    size = 3,
    fontface = "bold",
    na.rm = TRUE
  ) +
  
  ## rank: 1 at top, 10 at bottom
  scale_y_reverse(
    breaks = 1:10,
    labels = 1:10,
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  
  ## signatures on top
  scale_x_discrete(position = "top") +
  
  ## palettes
  scale_fill_manual(values = fill_pal, drop = FALSE, name = "FDR") +
  scale_color_identity() +
  
  ## we add tissue colors with a NEW scale by mapping tissue to color above,
  ## but color_identity() would override it — so we handle text color via I()
  NULL

## ---- IMPORTANT FIX:
## We must NOT use scale_color_identity() anymore because we need a color legend for tissues.
## So we rebuild the plot with explicit text color mapping via after_scale / I() approach.

p_heat_cat <- ggplot(hm_df, aes(x = signature, y = y_num, fill = fdr_bin)) +
  geom_tile(color = "white", linewidth = 0.6, na.rm = FALSE) +
  
  geom_rect(
    data = hm_df %>% filter(is_NR),
    inherit.aes = FALSE,
    aes(
      xmin = x_num - 0.5,
      xmax = x_num + 0.5,
      ymin = y_num + 0.5 - strip_h,
      ymax = y_num + 0.5
    ),
    fill = nr_purple,
    color = NA
  ) +
  
  ## dots: hollow background
  geom_point(
    data = dots_df,
    inherit.aes = FALSE,
    aes(x = x_dot, y = y_dot),
    shape  = 21,
    size   = dot_size,
    stroke = dot_stroke,
    fill   = "white",
    color  = "grey40",
    na.rm  = TRUE
  ) +
  ## dots: expressed (colored)
  geom_point(
    data = dots_df %>% filter(expressed),
    inherit.aes = FALSE,
    aes(x = x_dot, y = y_dot, color = tissue),
    shape  = 16,
    size   = dot_size,
    na.rm  = TRUE
  ) +
  
  ## TF labels with per-cell color:
  geom_text(
    aes(label = TF),
    color = hm_df$text_col,
    size = 3,
    fontface = "bold",
    na.rm = TRUE
  ) +
  
  scale_y_reverse(
    breaks = 1:10,
    labels = 1:10,
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_x_discrete(position = "top") +
  
  scale_fill_manual(values = fill_pal, drop = FALSE, name = "FDR") +
  scale_color_manual(
    name   = "TF expression (GTEx)",
    values = dot_pal,
    breaks = c("brain", "lung", "blood"),
    labels = c("Brain", "Lung", "Blood")
  ) +
  
  labs(x = NULL, y = "rank") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid  = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0),
    axis.text.y = element_text(face = "bold"),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, size = 1.2)
  )

p_heat_cat

