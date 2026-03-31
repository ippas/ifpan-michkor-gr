## ============================================================
## TF heatmap (DISCRETE FDR BINS) + nuclear receptor marker bar
## - if TF is in nuclearReceptors_GO0004879 -> add a purple strip at the bottom of the tile
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
## 0) INPUT: your nuclear receptor list
## ----------------------------
# Example: replace with YOUR vector (can be symbols like "NR3C1", "PPARG", "LXR", etc.)
# nuclearReceptors_GO0004879 <- c("NR3C1","NR3C2","PPARG","RXRA","RXRB","LXRA","LXRB")
# If you already have it defined above, keep it and do NOT redefine.

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
## 2) FDR bins (as requested)
## (0.01, 1e-5] ; (1e-5, 1e-10] ; (1e-10, 1e-20] ; (1e-20, 1e-30] ; ≤1e-30
## plus ">0.01"
## ----------------------------
fdr_breaks <- c(0, 1e-30, 1e-20, 1e-10, 1e-5, 1e-2, Inf)

fdr_labels <- c(
  "≤1e-30",
  "(1e-30, 1e-20]",
  "(1e-20, 1e-10]",
  "(1e-10, 1e-5]",
  "(1e-5, 0.01]",
  ">0.01"
)

fill_pal <- c(
  ">0.01"          = "white",
  "(1e-5, 0.01]"   = "#fae8e6",
  "(1e-10, 1e-5]"  = "#eba8a1",
  "(1e-20, 1e-10]" = "#d46b66",
  "(1e-30, 1e-20]" = "#b01728",
  "≤1e-30"         = "#5a0011"
)

text_pal <- c(
  ">0.01"          = "black",
  "(1e-5, 0.01]"   = "black",
  "(1e-10, 1e-5]"  = "black",
  "(1e-20, 1e-10]" = "white",
  "(1e-30, 1e-20]" = "white",
  "≤1e-30"         = "white"
)

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

## ----------------------------
## 4) Heatmap dataframe (complete grid + bins + coordinates)
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

    ## numeric coords for drawing rectangles
    x_num     = as.integer(signature),
    y_num     = as.integer(rank)
  )

## ----------------------------
## 5) Plot (discrete heatmap + purple bottom strip if NR)
## ----------------------------
nr_purple <- "#7B2CBF"   # violet/purple strip color
strip_h   <- 0.10        # strip thickness (0.10 = thin)

p_heat_cat <- ggplot(hm_df, aes(x = signature, y = y_num, fill = fdr_bin)) +
  ## base tiles
  geom_tile(color = "white", linewidth = 0.6, na.rm = FALSE) +

  ## nuclear receptor bottom strip (drawn inside the same tile)
  geom_rect(
    data = hm_df %>% filter(is_NR),
    inherit.aes = FALSE,
    aes(
      xmin = x_num - 0.5,
      xmax = x_num + 0.5,
      ymin = y_num - 0.5,
      ymax = y_num - 0.5 + strip_h
    ),
    fill = nr_purple,
    color = NA
  ) +

  ## TF label (bold) with per-bin text color
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

  ## discrete fill + text color
  scale_fill_manual(values = fill_pal, drop = FALSE, name = "FDR") +
  scale_color_identity() +

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
