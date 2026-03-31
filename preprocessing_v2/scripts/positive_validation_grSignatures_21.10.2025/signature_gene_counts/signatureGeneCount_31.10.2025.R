# ============================================================
# 📦 Pakiety
# ============================================================
library(ggplot2)
library(dplyr)
library(stringr)
library(patchwork)

# ============================================================
# 🧩 Funkcja uniwersalna: barplot liczby genów w sygnaturach
# ============================================================
plot_signature_gene_counts <- function(signature_keys, signature_labels, title) {
  df <- tibble(
    signature_key = signature_keys,
    signature_label = signature_labels,
    n_genes = sapply(flat_allGrSignatures_31.10.2025[signature_keys], length)
  ) %>%
    mutate(
      direction = if_else(str_detect(signature_label, "Up"), "Up", "Down"),
      signature_name = str_remove(signature_label, "(Up|Down)$"),
      signature_name = case_when(
        str_detect(signature_name, "globalGr") ~ "GlobalSignatures",
        TRUE ~ signature_name
      ),
      signature_name = factor(
        signature_name,
        levels = c("GlobalSignatures", "NeuralCells", "BloodCells", "LungCells")
      ),
      direction = factor(direction, levels = c("Up", "Down"))
    )
  
  df_sum <- df %>%
    group_by(signature_name) %>%
    summarise(total_genes = sum(n_genes),
              label_y = max(n_genes) + 10,
              .groups = "drop")
  
  ggplot(df, aes(x = signature_name, y = n_genes, fill = direction)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black") +
    geom_text(
      data = df_sum,
      aes(x = signature_name, y = label_y, label = total_genes),
      inherit.aes = FALSE,
      size = 4,
      fontface = "bold",
      color = "black"
    ) +
    scale_fill_manual(values = c("Up" = "#8b0000", "Down" = "#08306b")) +
    labs(
      title = title,
      x = NULL,
      y = "Number of genes",
      fill = "Regulation"
    ) +
    theme_classic(base_size = 13) +
    theme(
      text = element_text(color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black", face = "bold"),
      axis.text.x = element_text(angle = 35, hjust = 1, color = "black"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold")
    )
}

# ============================================================
# 📊 1️⃣ Full (oryginalne sygnatury)
# ============================================================
keys_full <- c(
  "global_GR_genes_globalUp5TissuesDerivedCells",
  "global_GR_genes_globalDown5TissuesDerivedCells",
  "NeuralCellsUp", "NeuralCellsDown",
  "BloodCellsUp", "BloodCellsDown",
  "LungCellsUp", "LungCellsDown"
)
labels_full <- c(
  "globalGrUp_5TissuesDerived", "globalGrDown_5TissuesDerived",
  "NeuralCellsUp", "NeuralCellsDown",
  "BloodCellsUp", "BloodCellsDown",
  "LungCellsUp", "LungCellsDown"
)

p_full <- plot_signature_gene_counts(keys_full, labels_full, "Full signatures")

# ============================================================
# 📊 2️⃣ Minus cluster BMC
# ============================================================
keys_cluster <- c(
  "minusClustersKPO_NeuralCellsUp",
  "minusClusterD_NeuralCellsDown",
  "minusClustersKPO_BloodCellsUp",
  "minusClusterD_BloodCellsDown",
  "minusClustersKPO_LungCellsUp",
  "minusClusterD_LungCellsDown"
)
labels_cluster <- c(
  "NeuralCellsUp", "NeuralCellsDown",
  "BloodCellsUp", "BloodCellsDown",
  "LungCellsUp", "LungCellsDown"
)

p_cluster <- plot_signature_gene_counts(keys_cluster, labels_cluster, "Minus cluster BMC")

# ============================================================
# 📊 3️⃣ Minus global (6 tissues)
# ============================================================
keys_6 <- c(
  "global_GR_genes_globalUp6TissuesDerivedCells",
  "global_GR_genes_globalDown6TissuesDerivedCells",
  "minusGlobalUp6TissuesDerivedCells_NeuralCellsUp",
  "minusGlobalDown6TissuesDerivedCells_NeuralCellsDown",
  "minusGlobalUp6TissuesDerivedCells_BloodCellsUp",
  "minusGlobalDown6TissuesDerivedCells_BloodCellsDown",
  "minusGlobalUp6TissuesDerivedCells_LungCellsUp",
  "minusGlobalDown6TissuesDerivedCells_LungCellsDown"
)
labels_6 <- c(
  "globalGrUp_6TissuesDerived", "globalGrDown_6TissuesDerived",
  "NeuralCellsUp", "NeuralCellsDown",
  "BloodCellsUp", "BloodCellsDown",
  "LungCellsUp", "LungCellsDown"
)

p6 <- plot_signature_gene_counts(keys_6, labels_6, "Minus global (6 tissues)")

# ============================================================
# 📊 4️⃣ Minus global (UP+DOWN, 5 tissues) — nowe listy
# ============================================================
keys_5 <- c(
  "global_GR_genes_globalUp5TissuesDerivedCells",
  "global_GR_genes_globalDown5TissuesDerivedCells",
  "minusGlobalUpDown5TissuesDerivedCells_NeuralCellsUp",
  "minusGlobalUpDown5TissuesDerivedCells_NeuralCellsDown",
  "minusGlobalUpDown5TissuesDerivedCells_BloodCellsUp",
  "minusGlobalUpDown5TissuesDerivedCells_BloodCellsDown",
  "minusGlobalUpDown5TissuesDerivedCells_LungCellsUp",
  "minusGlobalUpDown5TissuesDerivedCells_LungCellsDown"
)
labels_5 <- c(
  "globalGrUp_5TissuesDerived", "globalGrDown_5TissuesDerived",
  "NeuralCellsUp", "NeuralCellsDown",
  "BloodCellsUp", "BloodCellsDown",
  "LungCellsUp", "LungCellsDown"
)

p5 <- plot_signature_gene_counts(keys_5, labels_5, "Minus global (UP+DOWN, 5 tissues)")

# ============================================================
# 🧩 Połączenie 4 wykresów w układzie 2×2
# ============================================================
combined_plot <- ((p_full | p_cluster) / (p6 | p5)) +
  plot_annotation(title = "Comparison of GR-dependent signatures – gene counts (31.10.2025)") &
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

# ============================================================
# 💾 Zapis SVG
# ============================================================
output_dir <- "results_v2/positive_validation_grSignatures_31.10.2025/plots_signature_gene_counts"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

output_file <- file.path(output_dir, "GR_signatures_geneCounts_allVariants_31.10.2025.svg")

svg(output_file, width = 16, height = 12)
print(combined_plot)
dev.off()

message("✅ Zapisano plik: ", output_file)
