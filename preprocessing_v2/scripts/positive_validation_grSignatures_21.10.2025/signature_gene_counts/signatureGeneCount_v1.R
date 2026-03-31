library(ggplot2)
library(dplyr)
library(stringr)

# 🔹 Dane wejściowe
signature_keys <- c(
  "global_GR_genes_globalUp5TissuesDerivedCells",
  "global_GR_genes_globalDown5TissuesDerivedCells",
  "minusGlobalUp5TissuesDerivedCells_NeuralCellsUp",
  "minusGlobalDown5TissuesDerivedCells_NeuralCellsDown",
  "minusGlobalUp5TissuesDerivedCells_BloodCellsUp",
  "minusGlobalDown5TissuesDerivedCells_BloodCellsDown",
  "minusGlobalUp5TissuesDerivedCells_LungCellsUp",
  "minusGlobalDown5TissuesDerivedCells_LungCellsDown"
)

signature_labels <- c(
  "globalGrUp_5TissuesDerived",
  "globalGrDown_5TissuesDerived",
  "NeuralCellsUp",
  "NeuralCellsDown",
  "BloodCellsUp",
  "BloodCellsDown",
  "LungCellsUp",
  "LungCellsDown"
)

# 🔹 Tworzenie dataframe
gene_counts_df <- tibble(
  signature_key = signature_keys,
  signature_label = signature_labels,
  n_genes = sapply(flat_allGrSignatures_17.10.2025[signature_keys], length)
) %>%
  mutate(
    direction = if_else(str_detect(signature_label, "Up"), "Up-regulated", "Down-regulated"),
    signature_name = str_remove(signature_label, "(Up|Down)$"),
    signature_name = case_when(
      str_detect(signature_name, "globalGr") ~ "Global signatures",
      TRUE ~ signature_name
    ),
    signature_name = factor(signature_name, levels = c("Global signatures", "NeuralCells", "BloodCells", "LungCells")),
    direction = factor(direction, levels = c("Up-regulated", "Down-regulated"))
  )

# 🔹 Oblicz sumy i pozycję etykiet
gene_counts_sum <- gene_counts_df %>%
  group_by(signature_name) %>%
  summarise(
    total_genes = sum(n_genes),
    label_y = max(n_genes) + 20, # większy margines
    .groups = "drop"
  )

# 📁 Ścieżka zapisu
output_path <- "/home/mateusz/projects/ifpan-michkor-gr/results_v2/positive_validation_grSignatures_21.10.2025/signature_gene_counts/grSignaturesCounts_global5TissuesDerived_tissuesMinusGlobal_24.10.2025.svg"

# 🖼️ Otwarcie urządzenia graficznego SVG
svg(
  filename = output_path,
  width = 9,    # szerokość w calach
  height = 9,   # wysokość w calach
  pointsize = 12,  # rozmiar bazowy czcionki
  bg = "white"     # tło wykresu
)

# 🔹 Rysowanie wykresu
ggplot(gene_counts_df, aes(x = signature_name, y = n_genes, fill = direction)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black") +
  geom_text(
    data = gene_counts_sum,
    aes(x = signature_name, y = label_y, label = total_genes),
    inherit.aes = FALSE,
    size = 6,
    fontface = "bold",
    color = "black"
  ) +
  scale_fill_manual(values = c("Up-regulated" = "firebrick", "Down-regulated" = "darkblue")) +
  labs(
    x = "GR-dependent gene signatures",
    y = "Number of genes",
    fill = "Regulation direction"
  ) +
  theme_classic(base_size = 18) +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(color = "black", size = 16, face = "bold"),
    axis.title = element_text(color = "black", face = "bold", size = 18),
    legend.text = element_text(color = "black", size = 16),
    legend.title = element_text(color = "black", face = "bold", size = 17),
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 16, face = "bold"),
    legend.position = "bottom",
    legend.key.size = unit(1.2, "cm"),
    plot.margin = margin(10, 10, 10, 10)
  )

# ✅ Zamknięcie urządzenia graficznego
dev.off()

message("✅ Plik SVG zapisany do: ", output_path)
           
