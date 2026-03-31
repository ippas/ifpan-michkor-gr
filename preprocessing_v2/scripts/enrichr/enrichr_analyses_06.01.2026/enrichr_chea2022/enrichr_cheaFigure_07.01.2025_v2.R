## ============================================================
## TF heatmap (DISCRETE FDR BINS)
## - build top10 per signature (best term per TF = min FDR)
## - x = signature (ordered, shown on top)
## - y = rank (1 at top -> 10 at bottom)
## - fill = FDR category bins
## - text color (black/white) chosen per bin for readability
## - black frame around heatmap
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
## SIGNATURE ORDER
## ----------------------------
signature_order <- c(
  "global_up", "global_down",
  "brain_up",  "brain_down",
  "blood_up",  "blood_down",
  "lung_up",   "lung_down"
)

## ----------------------------
## FDR BINS (as requested)
## 0.01 -> 1e-10 -> 1e-20 -> 1e-30 -> 1e-40 -> smaller
## plus ">0.01" as last bin
## ----------------------------
fdr_breaks <- c(0, 1e-40, 1e-30, 1e-20, 1e-10, 1e-2, Inf)
fdr_labels <- c(
  "≤1e-40",
  "(1e-40, 1e-30]",
  "(1e-30, 1e-20]",
  "(1e-20, 1e-10]",
  "(1e-10, 0.01]",
  ">0.01"
)

## ----------------------------
## PALETTE: fill colors per bin
## (darker = stronger / smaller FDR)
## ----------------------------
fill_pal <- c(
  ">0.01"          = "white",
  "(1e-10, 0.01]"  = "#fae8e6",
  "(1e-20, 1e-10]" = "#eba8a1",
  "(1e-30, 1e-20]" = "#d46b66",
  "(1e-40, 1e-30]" = "#b01728",
  "≤1e-40"         = "#5a0011"
)

## ----------------------------
## TEXT COLOR per bin (readable)
## ----------------------------
text_pal <- c(
  ">0.01"          = "black",
  "(1e-10, 0.01]"  = "black",
  "(1e-20, 1e-10]" = "black",
  "(1e-30, 1e-20]" = "white",
  "(1e-40, 1e-30]" = "white",
  "≤1e-40"         = "white"
)

## ============================================================
## 1) BUILD ONE TABLE: top10 per signature
## ============================================================
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

## ============================================================
## UPDATED FDR BINS (as requested):
## (0.01, 1e-5] ; (1e-5, 1e-10] ; (1e-10, 1e-20] ; (1e-20, 1e-30] ; ≤1e-30
## plus ">0.01" as last bin
## ============================================================

nuclearReceptors_GO0004879 <- c(
  "ABHD2", "AHR", "AHRR", "AR", "ARNT", "ESR1", "ESR2", "ESRRA", "ESRRB", 
  "ESRRG", "GPER1", "HNF4A", "HNF4G", "NKX3-1", "NR1D1", "NR1D2", "NR1H2", 
  "NR1H3", "NR1H4", "NR1H5", "NR1I2", "NR1I3", "NR2C1", "NR2C2", "NR2E1", 
  "NR2E3", "NR2F1", "NR2F2", "NR2F6", "NR3C1", "NR3C2", "NR4A1", "NR4A2", 
  "NR4A3", "NR5A1", "NR5A2", "NR6A1", "OR51E2", "PAQR7", "PAQR8", "PDE3A", 
  "PGR", "POU5F1", "PPARA", "PPARD", "PPARG", "RARA", "RARB", "RARG", "RORA", 
  "RORB", "RORC", "RXRA", "RXRB", "RXRG", "SREBF1", "STAT3", "THRA", "THRB", 
  "VDR", "LXR"
)

fdr_breaks <- c(0, 1e-30, 1e-20, 1e-10, 1e-6, 1e-2, Inf)

fdr_labels <- c(
  "≤1e-30",
  "(1e-30, 1e-20]",
  "(1e-20, 1e-10]",
  "(1e-10, 1e-6]",
  "(1e-5, 0.01]",
  ">0.01"
)

fill_pal <- c(
  ">0.01"          = "white",
  "(1e-5, 0.01]"   = "#fae8e6",
  "(1e-10, 1e-6]"  = "#eba8a1",
  "(1e-20, 1e-10]" = "#d46b66",
  "(1e-30, 1e-20]" = "#b01728",
  "≤1e-30"         = "#5a0011"
)

## text color (black on light, white on dark)
text_pal <- c(
  ">0.01"          = "black",
  "(1e-5, 0.01]"   = "black",
  "(1e-10, 1e-6]"  = "black",
  "(1e-20, 1e-10]" = "black",
  "(1e-30, 1e-20]" = "white",
  "≤1e-30"         = "white"
)

## ============================================================
## Rebuild hm_df (only this part needs rerun from your pipeline)
## ============================================================
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
    text_col  = text_pal[as.character(fdr_bin)]
  )

p_heat_cat <- ggplot(hm_df, aes(x = signature, y = rank, fill = fdr_bin)) +
  geom_tile(color = "white", linewidth = 0.6, na.rm = FALSE) +
  geom_text(aes(label = TF, color = text_col), size = 12/2.83, fontface = "bold", na.rm = TRUE) +
  scale_y_discrete(limits = rev(levels(hm_df$rank))) +
  scale_x_discrete(position = "top") +
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
  # geom_text(
  #   aes(label = TF, color = text_col),
  #   size = 4,
  #   fontface = "bold",
  #   na.rm = TRUE
  # )

p_heat_cat
