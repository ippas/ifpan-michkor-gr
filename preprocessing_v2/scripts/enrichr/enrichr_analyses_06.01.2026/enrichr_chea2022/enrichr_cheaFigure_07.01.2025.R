## ============================================================
## TF heatmap: top10 per signature (best term per TF by min FDR)
## X = signature (ordered, shown on top)
## Y = rank (1 at top -> 10 at bottom)
## fill = -log10(FDR), white -> burgundy (darker = stronger)
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
})

signature_order <- c(
  "global_up", "global_down",
  "brain_up",  "brain_down",
  "blood_up",  "blood_down",
  "lung_up",   "lung_down"
)

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

hm_df <- top10_per_signature %>%
  mutate(
    signature = factor(signature, levels = signature_order),
    rank      = factor(rank, levels = 1:10, ordered = TRUE),
    neglog10  = -log10(pmax(FDR, 1e-300))
  ) %>%
  complete(signature, rank) %>%
  mutate(rank = factor(rank, levels = 1:10, ordered = TRUE))

p_heat <- ggplot(hm_df, aes(x = signature, y = rank, fill = neglog10)) +
  geom_tile(color = "white", linewidth = 0.6, na.rm = FALSE) +
  geom_text(aes(label = TF), size = 3, na.rm = TRUE) +
  scale_y_discrete(limits = rev(levels(hm_df$rank))) +  # 1 at top
  scale_fill_gradient(
    low  = "white",
    high = "#7a0019",
    name = expression(-log[10](FDR))
  ) +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = "rank") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid  = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0),
    axis.text.y = element_text(face = "bold"),
    legend.position = "right"
  )

p_heat +
  theme(
    panel.background = element_rect(colour = "black", fill = NA, size = 2)
  )
  
