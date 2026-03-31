# ============================================================
# 📦 Pakiety
# ============================================================
library(ggplot2)
library(dplyr)
library(stringr)
library(patchwork)

# ============================================================
# 🧠 Funkcja rysująca barplot dla dowolnego zestawu sygnatur
# ============================================================
plot_signature_gene_counts <- function(signature_keys, signature_labels, title) {
  df <- tibble(
    signature_key   = signature_keys,
    signature_label = signature_labels,
    n_genes         = sapply(flat_allGrSignatures_18.11.2025[signature_keys], length)
  ) %>%
    mutate(
      direction = if_else(str_detect(signature_label, "Up"), "Up", "Down"),
      signature_name = str_remove(signature_label, "(Up|Down)$"),
      signature_name = case_when(
        str_detect(signature_name, "globalGr") ~ "GlobalSignatures",
        TRUE ~ signature_name
      ),
      signature_name = factor(signature_name,
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
      aes(label = n_genes),
      position = position_dodge(width = 0.8),
      vjust = -0.5,
      size = 4,
      fontface = "bold"
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
      axis.text.x = element_text(angle = 35, hjust = 1, color = "black"),
      legend.position = "none"
    )
}

# ============================================================
# 1️⃣ RAW signatures (FULL)
# ============================================================
keys_full <- c(
  "global_GR_genes_globalUp",
  "global_GR_genes_globalDown",
  "NeuralCellsUp", "NeuralCellsDown",
  "BloodCellsUp", "BloodCellsDown",
  "LungCellsUp", "LungCellsDown"
)
labels_full <- c(
  "globalGrUp_RAW", "globalGrDown_RAW",
  "NeuralCellsUp", "NeuralCellsDown",
  "BloodCellsUp", "BloodCellsDown",
  "LungCellsUp", "LungCellsDown"
)

p_full <- plot_signature_gene_counts(keys_full, labels_full, "Full signatures")

# ============================================================
# 2️⃣ Minus GLOBAL (4 tissues)
# ============================================================
keys_4 <- c(
  "global_GR_genes_globalUp4TissuesDerivedCells",
  "global_GR_genes_globalDown4TissuesDerivedCells",
  "minusGlobalUpDown4TissuesDerivedCells_NeuralCellsUp",
  "minusGlobalUpDown4TissuesDerivedCells_NeuralCellsDown",
  "minusGlobalUpDown4TissuesDerivedCells_BloodCellsUp",
  "minusGlobalUpDown4TissuesDerivedCells_BloodCellsDown",
  "minusGlobalUpDown4TissuesDerivedCells_LungCellsUp",
  "minusGlobalUpDown4TissuesDerivedCells_LungCellsDown"
)
labels_4 <- c(
  "globalGrUp_4TissuesDerived", "globalGrDown_4TissuesDerived",
  "NeuralCellsUp", "NeuralCellsDown",
  "BloodCellsUp", "BloodCellsDown",
  "LungCellsUp", "LungCellsDown"
)

p4 <- plot_signature_gene_counts(keys_4, labels_4, "Minus global (4 tissues)")

# ============================================================
# 3️⃣ Minus GLOBAL (5 tissues)
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

p5 <- plot_signature_gene_counts(keys_5, labels_5, "Minus global (5 tissues)")

# ============================================================
# 4️⃣ Minus GLOBAL (6 tissues)
# ============================================================
keys_6 <- c(
  "global_GR_genes_globalUp6TissuesDerivedCells",
  "global_GR_genes_globalDown6TissuesDerivedCells",
  "minusGlobalUpDown6TissuesDerivedCells_NeuralCellsUp",
  "minusGlobalUpDown6TissuesDerivedCells_NeuralCellsDown",
  "minusGlobalUpDown6TissuesDerivedCells_BloodCellsUp",
  "minusGlobalUpDown6TissuesDerivedCells_BloodCellsDown",
  "minusGlobalUpDown6TissuesDerivedCells_LungCellsUp",
  "minusGlobalUpDown6TissuesDerivedCells_LungCellsDown"
)
labels_6 <- c(
  "globalGrUp_6TissuesDerived", "globalGrDown_6TissuesDerived",
  "NeuralCellsUp", "NeuralCellsDown",
  "BloodCellsUp", "BloodCellsDown",
  "LungCellsUp", "LungCellsDown"
)

p6 <- plot_signature_gene_counts(keys_6, labels_6, "Minus global (6 tissues)")

# ============================================================
# 🧩 COMBINE (2 × 2)
# ============================================================
final_plot <- ((p_full | p4) / (p5 | p6)) +
  plot_annotation(title = "Comparison of GR-dependent signatures – gene counts (18.11.2025)") &
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))

# ============================================================
# 💾 SAVE SVG
# ============================================================
output_dir <- "/home/mateusz/projects/ifpan-michkor-gr/results_v2/positive_validation_grSignatures_18.11.2025/plots_signature_gene_counts"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

output_file <- file.path(output_dir, "GR_signatures_geneCounts_full_minus4_minus5_minus6_18.11.2025.svg")

svg(output_file, width = 16, height = 12)
print(final_plot)
dev.off()

message("🔥 Zapisano plik: ", output_file)

