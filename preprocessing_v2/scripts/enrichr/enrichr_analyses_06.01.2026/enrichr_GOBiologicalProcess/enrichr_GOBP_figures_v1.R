## ============================================================
## GO BP barplot (UP) – group by tissue, no facet
## INPUT:
##   res$enrichr$overlap_only  -> named list of data.frames (per signature)
##   each df contains at least: Term, FDR
## OUTPUT:
##   ggplot object p_up
##   (optionally) saves SVG/PDF for Inkscape
## ============================================================
res <- run_enrichr_multi(
  gene_lists = flat_allGrSignatures_31.10.2025[sig_names],
  database = "GO_Biological_Process_2025",
  min_overlap_genes = 3,
  fdr_threshold = 0.05,
  names_mapper = names_mapper,
  xlsx_file = "results_v2/overlap/enrichr/enrichr_GOBiologicalProcess_GRsignatures_06.01.2026.xlsx"
)


tissue_colors <- c(
  global = "#9ECAE1",  # pastelowy niebieski
  brain  = "#FDD0A2",  # jasny pomarańcz
  blood  = "#FCAEAE",  # jasna czerwień
  lung   = "#A1D99B"   # jasna zieleń
)

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(stringr)
  library(forcats)
  library(ggplot2)
})

## -----------------------------
## SETTINGS
## -----------------------------
signatures_up <- c("global_up", "brain_up", "blood_up", "lung_up")
tissue_order  <- c("global", "brain", "blood", "lung")

# tissue_colors <- c(
#   global = "#4D4D4D",
#   brain  = "#6A51A3",
#   blood  = "#B2182B",
#   lung   = "#1B9E77"
# )

top_n_per_signature <- 10
fdr_line            <- 0.05

## optional saving (set to NULL to skip)
svg_file <- NULL  # e.g. "GO_BP_UP_top10_grouped.svg"
pdf_file <- NULL  # e.g. "GO_BP_UP_top10_grouped.pdf"

## -----------------------------
## LABEL FORMATTER (GO:ID + spaces + TERM)
## -----------------------------
format_go_label_id_first <- function(x, n_spaces = 5) {
  id   <- sub("^.*\\((GO:[0-9]+)\\)$", "\\1", x)      # extract "GO:xxxxxxx"
  name <- sub("\\s*\\(GO:[0-9]+\\)$", "", x)          # remove "(GO:xxxxxxx)"
  paste0(id, strrep(" ", n_spaces), name)             # GO:id + 5 spaces + term
}

## -----------------------------
## 1) TAKE TOP N FROM EACH SIGNATURE (UP)
## -----------------------------
up_list <- res$enrichr$overlap_only %>%
  lapply(head, top_n_per_signature) %>%
  .[signatures_up]

stopifnot(is.list(up_list), length(up_list) == length(signatures_up))

## -----------------------------
## 2) FLATTEN TO ONE DF + ANNOTATE
## -----------------------------
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
    tissue = factor(tissue, levels = tissue_order)
  )

## validate required columns
required_cols <- c("Term", "FDR", "n_genes", "signature", "tissue")
missing_cols  <- setdiff(required_cols, colnames(up_df))
if (length(missing_cols) > 0) {
  stop("Missing required columns in up_df: ", paste(missing_cols, collapse = ", "))
}

## score for x-axis
up_df <- up_df %>%
  mutate(
    FDR     = as.numeric(FDR),
    score   = -log10(FDR),
    n_genes = as.integer(n_genes)
  )

## -----------------------------
## 3) ORDER Y-AXIS: blocks by tissue, within each by score desc
## -----------------------------
up_df2 <- up_df %>%
  mutate(term_label = Term) %>%
  arrange(tissue, desc(score), term_label) %>%
  mutate(
    term_id = paste0(term_label, "___", tissue),
    term_id = factor(term_id, levels = rev(unique(term_id)))
  )

## -----------------------------
## 4) PLOT
## -----------------------------
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
    labels = function(x) format_go_label_id_first(sub("___.*$", "", x), n_spaces = 5)
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
    
    ## all text black
    text        = element_text(color = "black"),
    axis.text   = element_text(color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title  = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    
    axis.line   = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill  = NA,
      size  = 0.9
    )
  ) +
  expand_limits(x = max(up_df2$score, na.rm = TRUE) * 1.15)

print(p_up)


## ============================================================
## GO BP barplot (DOWN) – grouped by tissue, mirrored X-axis
##  - Y labels on the RIGHT
##  - X axis mirrored (bars go to the LEFT)
## INPUT:
##   res$enrichr$overlap_only  -> named list of data.frames (per signature)
##   each df contains at least: Term, FDR, n_genes
## OUTPUT:
##   ggplot object p_down
##   (optionally) saves SVG/PDF for Inkscape
## ============================================================


suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(stringr)
  library(forcats)
  library(ggplot2)
})

## -----------------------------
## SETTINGS
## -----------------------------
signatures_down <- c("global_down", "brain_down", "blood_down", "lung_down")
tissue_order    <- c("global", "brain", "blood", "lung")

# tissue_colors <- c(
#   global = "#4D4D4D",
#   brain  = "#6A51A3",
#   blood  = "#B2182B",
#   lung   = "#1B9E77"
# )

top_n_per_signature <- 10
fdr_line            <- 0.05

## optional saving (set to NULL to skip)
svg_file <- NULL  # e.g. "GO_BP_DOWN_top10_grouped_mirrored.svg"
pdf_file <- NULL  # e.g. "GO_BP_DOWN_top10_grouped_mirrored.pdf"

## -----------------------------
## LABEL FORMATTER (TERM + spaces + GO:ID)
## -----------------------------
format_go_label <- function(x, n_spaces = 5) {
  name <- sub("\\s*\\(GO:[0-9]+\\)$", "", x)          # remove "(GO:xxxxxxx)"
  id   <- sub("^.*\\((GO:[0-9]+)\\)$", "\\1", x)      # extract "GO:xxxxxxx"
  paste0(name, strrep(" ", n_spaces), id)             # name + 5 spaces + GO:id
}

## -----------------------------
## 1) TAKE TOP N FROM EACH SIGNATURE (DOWN)
## -----------------------------
down_list <- res$enrichr$overlap_only %>%
  lapply(head, top_n_per_signature) %>%
  .[signatures_down]

stopifnot(is.list(down_list), length(down_list) == length(signatures_down))

## -----------------------------
## 2) FLATTEN TO ONE DF + ANNOTATE
## -----------------------------
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
    tissue = factor(tissue, levels = tissue_order)
  )

## validate required columns
required_cols <- c("Term", "FDR", "n_genes", "signature", "tissue")
missing_cols  <- setdiff(required_cols, colnames(down_df))
if (length(missing_cols) > 0) {
  stop("Missing required columns in down_df: ", paste(missing_cols, collapse = ", "))
}

## score for significance (positive), and mirrored x (negative)
down_df <- down_df %>%
  mutate(
    FDR   = as.numeric(FDR),
    score = -log10(FDR),
    x     = -score,                # <-- MIRROR: bars go left
    n_genes = as.integer(n_genes)
  )

## -----------------------------
## 3) ORDER Y-AXIS: blocks by tissue, within each by score desc
## -----------------------------
down_df2 <- down_df %>%
  mutate(term_label = Term) %>%
  arrange(tissue, desc(score), term_label) %>%
  mutate(
    term_id = paste0(term_label, "___", tissue),
    term_id = factor(term_id, levels = rev(unique(term_id)))
  )

## -----------------------------
## 4) PREPARE LIMITS / THRESHOLD
## -----------------------------
x_min  <- min(down_df2$x, na.rm = TRUE)
x_pad  <- abs(x_min) * 0.15
x_line <- -(-log10(fdr_line))   # == -log10(0.05) with minus sign -> negative value (~ -1.301)

## -----------------------------
## 5) PLOT
## -----------------------------
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
    hjust = 1.15,      # place just beyond (left of) the bar end
    size  = 3,
    color = "black"
  ) +
  scale_y_discrete(
    labels = function(x) format_go_label(sub("___.*$", "", x), n_spaces = 5)
  ) +
  scale_x_continuous(
    labels = function(x) abs(x)  # show positive values on axis ticks
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
    
    ## all text black
    text        = element_text(color = "black"),
    axis.text   = element_text(color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title  = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    
    ## single clean border
    axis.line   = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill  = NA,
      size  = 0.9
    )
  ) +
  expand_limits(x = x_min - x_pad)  # extra space on the LEFT so labels don't clip

print(p_down)


p_down <- p_down +
  labs(subtitle = "Downregulated") +
  theme(
    plot.subtitle = element_text(
      hjust = 0.5,
      size  = 12,
      face  = "bold",
      color = "black",
      margin = margin(t = 15)   # <<< TO JEST KLUCZ
    )
  )


p_up <- p_up +
  labs(subtitle = "Upregulated") +
  theme(
    plot.subtitle = element_text(
      hjust = 0.5,
      size  = 12,
      face  = "bold",
      color = "black",
      margin = margin(t = 15)   # <<< TO JEST KLUCZ
    )
  )

final_plot <- (p_down | p_up) +
  plot_layout(guides = "collect") &
  theme(
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.title     = element_blank(),
    legend.text      = element_text(color = "black")
  )

final_plot

# tissue_colors <- c(
#   global = "#8DA0CB",  # pastelowy niebiesko-lawendowy
#   brain  = "#FC8D62",  # pastelowa morela
#   blood  = "#E78AC3",  # pastelowy róż
#   lung   = "#66C2A5"   # pastelowa mięta
# )
