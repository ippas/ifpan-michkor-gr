## ============================================================
## TF heatmap (DISCRETE FDR BINS) + nuclear receptor bottom strip
## + GTEx expression dots
## - FDR bins include 1e-6
## - nuclear receptor list: nuclearReceptors_GO0004879
## - purple strip at the BOTTOM of each tile (with y reversed)
## - dots:
##   * global_* : 3 dots (Brain/Lung/Blood)
##   * brain_*  : 1 dot (Brain only)
##   * lung_*   : 1 dot (Lung only)
##   * blood_*  : 1 dot (Blood only)
##   * hollow if NOT expressed; colored if expressed (median_max > 1)
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
## 1) Signature order
## ----------------------------
signature_order <- c(
  "global_up", "global_down",
  "brain_up",  "brain_down",
  "blood_up",  "blood_down",
  "lung_up",   "lung_down"
)

## ----------------------------
## 2) FDR bins
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
## 3) Top10 per signature
## ----------------------------
top10_per_signature <- res$enrichr$overlap_only %>%
  imap_dfr(~ .x %>%
             mutate(
               TF  = str_extract(Term, "^[^ ]+"),
               FDR = as.numeric(FDR)
             ) %>%
             group_by(TF) %>%
             slice_min(FDR, n = 1, with_ties = FALSE) %>%
             ungroup() %>%
             arrange(FDR) %>%
             slice_head(n = 10) %>%
             mutate(rank = row_number(),
                    signature = .y) %>%
             select(signature, rank, TF, everything())
  ) %>%
  mutate(
    TF = case_when(
      TF == "LXR"      ~ "NR1H3",
      TF == "SA1"      ~ "STAG1",
      TF == "SMC1"     ~ "SMC1A",
      TF == "CJUN"     ~ "JUN",
      TF == "TCFCP2L1" ~ "TFCP2L1",
      TF %in% c("Nerf2", "NRF2") ~ "NFE2L2",
      TF == "RING1B"   ~ "RING1",
      TF == "AF4"      ~ "AFF4",
      TRUE ~ TF
    )
  )

## ----------------------------
## 4) GTEx download
## ----------------------------
TF_expression_GTEx <- multi_downloadMedianSubtissuesGTEx_tissueSummaryExpression(
  gene_symbol   = top10_per_signature$TF,
  verbose       = TRUE,
  show_progress = TRUE,
  keep_unmapped = FALSE
)

## ----------------------------
## 5) Heatmap df
## ----------------------------
hm_df <- top10_per_signature %>%
  mutate(
    signature = factor(signature, levels = signature_order),
    rank      = factor(rank, levels = 1:10, ordered = TRUE)
  ) %>%
  complete(signature, rank) %>%
  mutate(
    FDR_num = as.numeric(FDR),
    fdr_bin = cut(
      FDR_num,
      breaks = fdr_breaks,
      labels = fdr_labels,
      include.lowest = TRUE
    ),
    fdr_bin = factor(fdr_bin, levels = rev(fdr_labels)),
    text_col = text_pal[as.character(fdr_bin)],
    is_NR = !is.na(TF) & TF %in% nuclearReceptors_GO0004879,
    x_num = as.integer(signature),
    y_num = as.integer(rank)
  )

## ----------------------------
## 6) GTEx flags
## ----------------------------
TF_expr_flags <- TF_expression_GTEx %>%
  filter(tissue %in% c("Brain", "Lung", "Whole_Blood")) %>%
  mutate(
    tissue = recode(
      tissue,
      "Brain" = "brain",
      "Lung" = "lung",
      "Whole_Blood" = "blood"
    ),
    expressed = median_max > 1
  ) %>%
  group_by(query_gene_symbol, tissue) %>%
  summarise(expressed = any(expressed), .groups = "drop") %>%
  pivot_wider(names_from = tissue, values_from = expressed, values_fill = FALSE)

## ----------------------------
## 7) Dots df
## ----------------------------
dots_df <- hm_df %>%
  mutate(signature_chr = as.character(signature)) %>%
  left_join(TF_expr_flags, by = c("TF" = "query_gene_symbol")) %>%
  mutate(
    brain = ifelse(is.na(brain), FALSE, brain),
    lung  = ifelse(is.na(lung),  FALSE, lung),
    blood = ifelse(is.na(blood), FALSE, blood),
    sig_group = case_when(
      str_detect(signature_chr, "^global_") ~ "global",
      str_detect(signature_chr, "^brain_")  ~ "brain",
      str_detect(signature_chr, "^lung_")   ~ "lung",
      str_detect(signature_chr, "^blood_")  ~ "blood"
    )
  ) %>%
  pivot_longer(c(brain, lung, blood),
               names_to = "tissue",
               values_to = "expressed") %>%
  mutate(
    x_dot = x_num + case_when(
      sig_group == "global" & tissue == "brain" ~ -0.23,
      sig_group == "global" & tissue == "lung"  ~  0.00,
      sig_group == "global" & tissue == "blood" ~  0.23,
      TRUE ~ 0
    ),
    y_dot = y_num - 0.30
  ) %>%
  filter(sig_group == "global" | tissue == sig_group)

## ----------------------------
## 8) Plot
## ----------------------------
tile_width  <- 0.97
tile_height <- 0.97

p_heat_cat <- ggplot(hm_df, aes(signature, y_num, fill = fdr_bin)) +
  geom_tile(
    width = tile_width,
    height = tile_height,
    color = "white",
    linewidth = 0.6
  ) +
  geom_rect(
    data = hm_df %>% filter(is_NR),
    inherit.aes = FALSE,
    aes(
      xmin = x_num - tile_width/2,
      xmax = x_num + tile_width/2,
      ymin = y_num + tile_height/2 - 0.10,
      ymax = y_num + tile_height/2
    ),
    fill = "#B388FF"
  ) +
  geom_point(
    data = dots_df,
    aes(x_dot, y_dot),
    inherit.aes = FALSE,
    shape = 21,
    size = 3,
    fill = "white",
    color = "grey40",
    stroke = 0.8
  ) +
  geom_point(
    data = dots_df %>% filter(expressed),
    aes(x_dot, y_dot, color = tissue),
    inherit.aes = FALSE,
    size = 3
  ) +
  geom_text(
    aes(label = TF),
    color = hm_df$text_col,
    size = 3,
    fontface = "bold"
  ) +
  scale_y_reverse(breaks = 1:10) +
  scale_x_discrete(position = "top") +
  scale_fill_manual(values = fill_pal, name = "FDR") +
  scale_color_manual(
    values = c(brain = "#2E86FF", lung = "#00A676", blood = "#E31A1C"),
    name = "TF expression (GTEx)"
  ) +
  labs(x = NULL, y = "rank") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 0),
    axis.text.y = element_text(face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, size = 1.2)
  )

p_heat_cat
