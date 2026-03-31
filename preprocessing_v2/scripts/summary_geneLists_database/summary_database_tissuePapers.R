# =========================
# GR-dependent genes across tissues
# Barplots: UP (right) & DOWN (left, mirrored)
# Sortowanie: rosnąco po UP
# Wszystkie teksty: czarne
# Oś X: ticki 0, 2000, 4000
# =========================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(forcats)
  library(tibble)
})

# -------------------------
# 1) Dane wejściowe
# -------------------------
df <- tribble(
  ~tissue,             ~UP,  ~DOWN,
  "Neural",            3592,  3974,
  "Lung",              2511,  2306,
  "Blood",             2881,  2950,
  "Adipose",            392,   167,
  "Embryos",           1206,  1115,
  "Skeletal Muscle",    766,   399,
  "Adrenal gland",     1081,   964,
  "Bone",               333,   243,
  "Kidney",            1045,  1135,
  "Cartilage",          408,   632
)

# -------------------------
# 2) Sortowanie tkanek
# -------------------------
df <- df %>%
  mutate(tissue = fct_reorder(tissue, UP, .desc = FALSE))

# -------------------------
# 3) Wspólny motyw: wszystko na czarno
# -------------------------
theme_black_text <- theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(colour = "black", face = "bold"),
    axis.title = element_text(colour = "black"),
    axis.text  = element_text(colour = "black")
  )

# -------------------------
# 4) Barplot UP (w prawo)
# -------------------------
p_up <- ggplot(df, aes(x = tissue, y = UP)) +
  geom_col(fill = "firebrick", width = 0.7) +
  geom_text(aes(label = UP), hjust = -0.15, size = 4, colour = "black") +
  coord_flip() +
  scale_y_continuous(
    breaks = c(0, 2000, 4000),
    limits = c(0, 4000),
    expand = expansion(mult = c(0, 0.15))
  ) +
  labs(
    title = "GR-dependent genes across tissues (UP)",
    x = "Tissue (incl. in vitro cultures)",
    y = "Number of genes"
  ) +
  theme_black_text

# -------------------------
# 5) Barplot DOWN (lustrzane w lewo)
# -------------------------
p_down <- ggplot(df, aes(x = tissue, y = -DOWN)) +
  geom_col(fill = "navyblue", width = 0.7) +
  geom_text(aes(label = DOWN), hjust = 1.15, size = 4, colour = "black") +
  coord_flip() +
  scale_y_continuous(
    breaks = c(-4000, -2000, 0),
    labels = function(x) abs(x),
    limits = c(-4000, 0),
    expand = expansion(mult = c(0.15, 0))
  ) +
  labs(
    title = "GR-dependent genes across tissues (DOWN)",
    x = "Tissue (incl. in vitro cultures)",
    y = "Number of genes"
  ) +
  theme_black_text

# -------------------------
# 6) Wyświetlenie
# -------------------------
print(p_up)
print(p_down)

# -------------------------
# 7) (Opcjonalnie) Zapis do plików
# -------------------------
# ggsave("gr_genes_UP_barplot.png", p_up, width = 7, height = 4, dpi = 300)
# ggsave("gr_genes_DOWN_barplot.png", p_down, width = 7, height = 4, dpi = 300)

p_down + p_up -> p_combined


# -------------------------
# 6) Save to SVG (custom path)
# -------------------------
dev.off()
out_dir <- "/home/mateusz/projects/ifpan-michkor-gr/results_v2/grDatabase_preprocessing/supplementaryFigures"
out_svg <- file.path(out_dir, "gr_genes_up_down_mirrored.svg")

svg(filename = out_svg, width = 5, height = 5)  # inches
print(p_combined)
dev.off()


# =========================
# GR-dependent genes across tissues
# One barplot: N papers
# Pastel brown bars, reversed Y-axis
# =========================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tibble)
})

# -------------------------
# 1) Dane
# -------------------------
df_papers <- tribble(
  ~tissue,             ~n_papers,
  "Neural",            14,
  "Lung",               7,
  "Blood",              6,
  "Adipose",            2,
  "Embryos",            2,
  "Skeletal Muscle",    2,
  "Adrenal gland",      1,
  "Bone",               1,
  "Kidney",             1,
  "Cartilage",          1
) %>%
  mutate(
    tissue = factor(tissue, levels = rev(tissue))  # odwrócona kolejność osi Y
  )

# -------------------------
# 2) Wykres
# -------------------------
p_papers <- ggplot(df_papers, aes(x = tissue, y = n_papers)) +
  geom_col(width = 0.7, fill = "#703c1e") +   # pastelowy brąz
  geom_text(aes(label = n_papers), hjust = -0.15, size = 4, colour = "black") +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "GR-dependent genes across tissues",
    x = "Tissue (incl. in vitro cultures)",
    y = "N papers"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(colour = "black", face = "bold"),
    axis.title = element_text(colour = "black"),
    axis.text  = element_text(colour = "black")
  )

# -------------------------
# 3) Wyświetlenie
# -------------------------
print(p_papers)

# -------------------------
# 4) (Opcjonalnie) Zapis do SVG
# -------------------------
dev.off()
out_svg <- "/home/mateusz/projects/ifpan-michkor-gr/results_v2/grDatabase_preprocessing/supplementaryFigures/gr_tissue_paper_counts.svg"
svg(out_svg, width = 3, height = 5)
print(p_papers)
dev.off()
