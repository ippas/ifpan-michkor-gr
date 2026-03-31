## ============================================================
## GO BP barplots (UP + DOWN) – grouped by tissue
## - UP: bars to the RIGHT
## - DOWN: mirrored bars to the LEFT (axis labels shown as positive)
## - Y labels: ONLY term name (NO GO:ID)
## - shows n_genes on bars
## - combines into one figure with shared legend
##
## INPUT:
##   res$enrichr$overlap_only -> named list of data.frames per signature
##   required columns in each df: Term, FDR, n_genes
##
## OUTPUT:
##   p_down, p_up, final_plot
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(patchwork)  # for (p_down | p_up) + plot_layout()
})

## -----------------------------
## COLORS (tissue)
## -----------------------------
tissue_colors <- c(
  global = "#9ECAE1",  # pastel blue
  brain  = "#FDD0A2",  # light orange
  blood  = "#FCAEAE",  # light red
  lung   = "#A1D99B"   # light green
)

## -----------------------------
## SETTINGS
## -----------------------------
top_n_per_signature <- 10
fdr_line            <- 0.05

signatures_up   <- c("global_up", "brain_up", "blood_up", "lung_up")
signatures_down <- c("global_down", "brain_down", "blood_down", "lung_down")
tissue_order    <- c("global", "brain", "blood", "lung")

## optional saving (set NULL to skip)
svg_file <- NULL
pdf_file <- NULL

## -----------------------------
## LABEL FORMATTER (TERM ONLY)
## -----------------------------
term_only <- function(x) {
  sub("\\s*\\(GO:[0-9]+\\)$", "", x)   # remove trailing " (GO:xxxxxxx)"
}

## -----------------------------
## 1) UP: take top N per signature
## -----------------------------
up_list <- res$enrichr$overlap_only %>%
  lapply(head, top_n_per_signature) %>%
  .[signatures_up]

stopifnot(is.list(up_list), length(up_list) == length(signatures_up))

up_df <- up_list %>%
  bind_rows(.id = "signature") %>%
  mutate(
    direction = "up",
    tissue = case_when(
      str_detect(signature, "^global_") ~ "global",
      str_detect(signature, "^brain_")  ~ "brain",
      str_detect(signature, "^blood_")  ~ "blood",
      str_detect(signature, "^lung_")   ~ "lung",
      TRUE ~ "other"
    ),
    tissue = factor(tissue, levels = tissue_order),
    FDR     = as.numeric(FDR),
    score   = -log10(FDR),
    n_genes = as.integer(n_genes)
  )

required_cols <- c("Term", "FDR", "n_genes", "signature", "tissue")
missing_cols  <- setdiff(required_cols, colnames(up_df))
if (length(missing_cols) > 0) {
  stop("Missing required columns in UP df: ", paste(missing_cols, collapse = ", "))
}

## order y-axis: tissue blocks, within by score desc
up_df2 <- up_df %>%
  mutate(term_label = term_only(Term)) %>%
  arrange(tissue, desc(score), term_label) %>%
  mutate(
    term_id = paste0(term_label, "___", tissue),
    term_id = factor(term_id, levels = rev(unique(term_id)))
  )

p_up <- ggplot(up_df2, aes(x = score, y = term_id, fill = tissue)) +
  geom_col(width = 0.8) +
  geom_vline(
    xintercept = -log10(fdr_line),
    linetype   = "dashed",
    size       = 0.4,
    color      = "grey40"
  ) +
  geom_text(
    aes(label = n_genes),
    hjust = -0.15,
    size  = 3,
    color = "black"
  ) +
  scale_y_discrete(
    position = "right",
    labels = function(x) sub("___.*$", "", x)
  ) +
  scale_fill_manual(values = tissue_colors, drop = FALSE) +
  guides(fill = guide_legend(nrow = 1)) +
  labs(
    x = expression(-log[10]("FDR")),
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.title     = element_blank(),
    text        = element_text(color = "black"),
    axis.text   = element_text(color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title  = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    axis.line   = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.9)
  ) +
  expand_limits(x = max(up_df2$score, na.rm = TRUE) * 1.15)

p_up <- p_up +
  labs(subtitle = "Upregulated") +
  theme(
    plot.subtitle = element_text(
      hjust = 0.5, size = 12, face = "bold", color = "black",
      margin = margin(t = 15)
    )
  )

## -----------------------------
## 2) DOWN: take top N per signature
## -----------------------------
down_list <- res$enrichr$overlap_only %>%
  lapply(head, top_n_per_signature) %>%
  .[signatures_down]

stopifnot(is.list(down_list), length(down_list) == length(signatures_down))

down_df <- down_list %>%
  bind_rows(.id = "signature") %>%
  mutate(
    direction = "down",
    tissue = case_when(
      str_detect(signature, "^global_") ~ "global",
      str_detect(signature, "^brain_")  ~ "brain",
      str_detect(signature, "^blood_")  ~ "blood",
      str_detect(signature, "^lung_")   ~ "lung",
      TRUE ~ "other"
    ),
    tissue = factor(tissue, levels = tissue_order),
    FDR     = as.numeric(FDR),
    score   = -log10(FDR),
    x       = -score,                  # mirror
    n_genes = as.integer(n_genes)
  )

missing_cols  <- setdiff(required_cols, colnames(down_df))
if (length(missing_cols) > 0) {
  stop("Missing required columns in DOWN df: ", paste(missing_cols, collapse = ", "))
}

down_df2 <- down_df %>%
  mutate(term_label = term_only(Term)) %>%
  arrange(tissue, desc(score), term_label) %>%
  mutate(
    term_id = paste0(term_label, "___", tissue),
    term_id = factor(term_id, levels = rev(unique(term_id)))
  )

x_min  <- min(down_df2$x, na.rm = TRUE)
x_pad  <- abs(x_min) * 0.15
x_line <- -(-log10(fdr_line))  # negative threshold

p_down <- ggplot(down_df2, aes(x = x, y = term_id, fill = tissue)) +
  geom_col(width = 0.8) +
  geom_vline(
    xintercept = x_line,
    linetype   = "dashed",
    size       = 0.4,
    color      = "grey40"
  ) +
  geom_text(
    aes(label = n_genes),
    hjust = 1.15,
    size  = 3,
    color = "black"
  ) +
  scale_y_discrete(
    labels = function(x) sub("___.*$", "", x)
  ) +
  scale_x_continuous(
    labels = function(x) abs(x)  # show positive tick labels
  ) +
  scale_fill_manual(values = tissue_colors, drop = FALSE) +
  guides(fill = guide_legend(nrow = 1)) +
  labs(
    x = expression(-log[10]("FDR")),
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.title     = element_blank(),
    text        = element_text(color = "black"),
    axis.text   = element_text(color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title  = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    axis.line   = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.9)
  ) +
  expand_limits(x = x_min - x_pad)

p_down <- p_down +
  labs(subtitle = "Downregulated") +
  theme(
    plot.subtitle = element_text(
      hjust = 0.5, size = 12, face = "bold", color = "black",
      margin = margin(t = 15)
    )
  )

## -----------------------------
## 3) COMBINE (shared legend)
## -----------------------------
final_plot <- (p_down | p_up) +
  plot_layout(guides = "collect") &
  theme(
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.title     = element_blank(),
    legend.text      = element_text(color = "black")
  )

final_plot
