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

library(dplyr)
library(tidyr)
library(ggplot2)

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
    )
  )

df_long <- df_long %>%
  mutate(
    signature_name = factor(
      signature_name,
      levels = c(
        "Lung_plus", "Lung_minus",
        "Blood_plus", "Blood_minus",
        "Brain_plus", "Brain_minus",
        "Systemic_plus", "Systemic_minus"
      )
    )
  )



df_long <- df_long %>%
  mutate(
    disorder = factor(
      disorder,
      levels = c("MDD", "SCZ", "AD")
    )
  )



p_disgenet <- df_long %>%
  filter(database == "DisGeNET") %>%
  ggplot(aes(x = n_genes, y = signature_name)) +
  geom_col(width = 0.7, fill = "#4B6FA5") +
  geom_text(
    aes(label = n_genes),
    hjust = 1.1,   # przesunięcie w lewo (bo oś X odwrócona)
    color = "black",
    size = 3
  ) +
  facet_wrap(~ disorder, ncol = 1, nrow = 3) +
  scale_x_reverse(position = "top") +
  labs(
    title = "DisGeNET",
    x = "Number of overlapping genes",
    y = "GR-dependent gene signature"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold", color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(color = "black"),
    axis.title.x = element_text(color = "black"),
    axis.title.y = element_text(color = "black"),
    plot.title = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black")
  )

p_disgenet


p_gwas <- df_long %>%
  filter(database == "GWAS Catalog") %>%
  ggplot(aes(x = n_genes, y = signature_name)) +
  geom_col(width = 0.7, fill = "#B44A4A") +
  geom_text(
    aes(label = n_genes),
    hjust = -0.2,
    color = "black",
    size = 4
  ) +
  facet_wrap(~ disorder, ncol = 1, nrow = 3) +
  scale_x_continuous(
    limits = c(0, 20),
    position = "top",
    expand = expansion(mult = c(0.05, 0.15))
  ) +
  scale_y_discrete(position = "right") +
  labs(
    title = "GWAS Catalog",
    x = "Number of overlapping genes",
    y = "GR-dependent gene signature"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold", color = "black"),
    axis.text.y.right = element_text(size = 10, color = "black"),
    axis.text.x = element_text(color = "black"),
    axis.title.x = element_text(color = "black"),
    axis.title.y.right = element_text(color = "black"),
    plot.title = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black")
  )

p_gwas


p_disgenet + p_gwas

