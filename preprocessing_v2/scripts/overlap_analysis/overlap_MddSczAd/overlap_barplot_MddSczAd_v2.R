# =========================
# 0) LIBS
# =========================
library(dplyr)
library(tidyr)
library(ggplot2)

# =========================
# 1) INPUT DATA (df)
# =========================
df <- data.frame(
  signature_name = c(
    "Systemic_plus",
    "Systemic_minus",
    "Brain_plus",
    "Brain_minus",
    "Blood_plus",
    "Blood_minus",
    "Lung_plus",
    "Lung_minus"
  ),
  MDD_disgenet = c(12, 5, 6, 4, 1, 3, 5, 7),
  MDD_GWASCatalog = c(13, 8, 10, 3, 3, 6, 10, 15),
  SCZ_disgenet = c(21, 15, 16, 7, 6, 6, 13, 16),
  SCZ_GWASCatalog = c(15, 9, 15, 5, 7, 8, 12, 16),
  AD_disgenet = c(8, 3, 3, 4, 5, 3, 3, 3),
  AD_GWASCatalog = c(14, 7, 7, 2, 3, 0, 8, 7),
  stringsAsFactors = FALSE
)

# =========================
# 2) COLORS: 4 tissues (Twoje kolory)
# =========================
tissue_colors <- c(
  Systemic = "#DE77AE",
  Blood    = "#5AAE61",
  Lung     = "#4292C6",
  Brain    = "#F4A582"
)

# =========================
# 3) LONG FORMAT
# =========================
df_long <- df %>%
  pivot_longer(
    cols = -signature_name,
    names_to = "source",
    values_to = "n_genes"
  ) %>%
  mutate(
    disorder = case_when(
      grepl("^MDD", source) ~ "MDD",
      grepl("^SCZ", source) ~ "SCZ",
      grepl("^AD", source)  ~ "AD"
    ),
    database = case_when(
      grepl("disgenet", source) ~ "DisGeNET",
      grepl("GWASCatalog", source) ~ "GWAS Catalog"
    ),
    disorder = factor(disorder, levels = c("MDD", "SCZ", "AD")),
    signature_name = factor(
      signature_name,
      levels = c(
        "Lung_minus", "Lung_plus",
        "Blood_minus", "Blood_plus",
        "Brain_minus", "Brain_plus",
        "Systemic_minus", "Systemic_plus"
      )
    ),
    tissue = case_when(
      grepl("^Systemic", signature_name) ~ "Systemic",
      grepl("^Brain", signature_name)    ~ "Brain",
      grepl("^Blood", signature_name)    ~ "Blood",
      grepl("^Lung", signature_name)     ~ "Lung"
    )
  )

# =========================
# 4) PLOT: DisGeNET (mirror)
# =========================
p_disgenet <- df_long %>%
  filter(database == "DisGeNET") %>%
  ggplot(aes(x = n_genes, y = signature_name, fill = tissue)) +
  geom_col(width = 0.7) +
  geom_text(
    aes(label = n_genes),
    hjust = 1.1,
    color = "black",
    size = 4
  ) +
  facet_wrap(~ disorder, ncol = 1, nrow = 3) +
  scale_fill_manual(values = tissue_colors, name = "Tissue") +
  scale_x_reverse(
    limits = c(22, 0),
    # position = "top",
    expand = expansion(mult = c(0.05, 0.15))
  ) +
  labs(
    title = "DisGeNET",
    x = "Number of overlapping genes",
    y = "GR-dependent gene signature"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.text   = element_text(face = "bold", color = "black"),
    axis.text.y  = element_text(size = 10, color = "black"),
    axis.text.x  = element_text(color = "black"),
    axis.title.x = element_text(color = "black"),
    axis.title.y = element_text(color = "black"),
    plot.title   = element_text(color = "black"),
    legend.title = element_text(color = "black"),
    legend.text  = element_text(color = "black")
  )

# =========================
# 5) PLOT: GWAS Catalog (Y axis on right)
# =========================
p_gwas <- df_long %>%
  filter(database == "GWAS Catalog") %>%
  ggplot(aes(x = n_genes, y = signature_name, fill = tissue)) +
  geom_col(width = 0.7) +
  geom_text(
    aes(label = n_genes),
    hjust = -0.1,
    color = "black",
    size = 4
  ) +
  facet_wrap(~ disorder, ncol = 1, nrow = 3) +
  scale_fill_manual(values = tissue_colors, name = "Tissue") +
  scale_x_continuous(
    limits = c(0, 20),
    # position = "top",
    expand = expansion(mult = c(0.05, 0.15))
  ) +
  scale_y_discrete(position = "right") +
  labs(
    title = "GWAS Catalog",
    x = "Number of overlapping genes",
    y = "GR-dependent gene signature"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.text         = element_text(face = "bold", color = "black"),
    axis.text.y.right = element_text(size = 10, color = "black"),
    axis.text.x       = element_text(color = "black"),
    axis.title.x      = element_text(color = "black"),
    axis.title.y.right= element_text(color = "black"),
    plot.title        = element_text(color = "black"),
    legend.title     = element_text(color = "black"),
    legend.text      = element_text(color = "black")
  )

# =========================
# 6) PRINT
# =========================

final_plot <- (p_disgenet + p_gwas) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

dev.off()

# zapis do SVG (base R)
svg(
  filename = "/home/mateusz/projects/ifpan-michkor-gr/results_v2/overlap/overlap_MddSczAd/overlap_MDD_SCZ_AD_DisGeNET_GWAS.svg",
  width = 6,   # cale ~ 180 mm
  height = 10  # cale ~ 120 mm
)

print(final_plot)
dev.off()
