library(tidyverse)

df <- tibble::tibble(
  icd10_category = c("F0x","F1x","F2x","F3x","F4x","F8x"),
  mean_cs        = c(2.98, 2.35, 4.33, 4.85, 2.28, 2.02)
)

# style (ten sam co u Ciebie)
bar_fill  <- "#335C67"  # earth-tone brown
bar_edge  <- "black"
bar_width <- 0.7
bar_lwd   <- 1

theme_black_text <- theme(
  plot.title   = element_text(color = "black"),
  axis.title.x = element_text(color = "black"),
  axis.title.y = element_text(color = "black"),
  axis.text.x  = element_text(color = "black"),
  axis.text.y  = element_text(color = "black"),
  axis.ticks   = element_line(color = "black"),
  axis.line    = element_line(color = "black")
)

# jeden barplot: Y = icd10_category, X = mean_cs, oś X na górze
p_mean_cs <- ggplot(
  df %>%
    complete(icd10_category = paste0("F", 0:9, "x"), fill = list(mean_cs = 0)) %>%
    arrange(mean_cs) %>%  # low -> high (high na górze)
    mutate(icd10_category = factor(icd10_category, levels = icd10_category)),
  aes(x = mean_cs, y = icd10_category)
) +
  geom_col(width = bar_width, fill = bar_fill, color = bar_edge, linewidth = bar_lwd, alpha = 1) +
  scale_x_continuous(position = "top") +
  theme_classic() +
  theme_black_text

p_mean_cs

dev.off()
out_dir  <- "/home/mateusz/projects/ifpan-michkor-gr/results_v2/overlap/overlap_GWASCatalogDisGeNETgenebass/figures"
out_file <- paste0("meanCS_grSystemic_barplot_", "06.03.2026", ".svg")

svg(filename = file.path(out_dir, out_file), width = 5, height = 4)
print(p_mean_cs)
dev.off()

