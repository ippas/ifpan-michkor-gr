# ============================================================
# 📦 Wczytanie pakietów
# ============================================================
library(dplyr)
library(purrr)
library(ggplot2)
library(tidyr)

# ============================================================
# 🧠 1. Uniwersalna funkcja do analizy Enrichr
# ============================================================
get_enrichr_results_generic <- function(
    tissue,
    database = "ChEA_2022",
    prefix_up = "",
    prefix_down = ""
) {
  # Dynamiczne składanie nazw
  up_name <- paste0(prefix_up, tissue, "Up")
  down_name <- paste0(prefix_down, tissue, "Down")
  
  # Pobranie list genów
  up_list <- flat_allGrSignatures_17.10.2025[[up_name]]
  down_list <- flat_allGrSignatures_17.10.2025[[down_name]]
  
  # Ostrzeżenia, jeśli listy nie istnieją
  if (is.null(up_list)) warning(paste("⚠️ Nie znaleziono listy:", up_name))
  if (is.null(down_list)) warning(paste("⚠️ Nie znaleziono listy:", down_name))
  
  list(
    Up = run_enrichr(gene_list = up_list, database = database),
    Down = run_enrichr(gene_list = down_list, database = database)
  )
}


# ============================================================
# 🧩 3. Uniwersalna funkcja do tworzenia wykresu TF (np. NR3C1/NR3C2)
# ============================================================
plot_enrichr_tf_results <- function(
    enrichr_results,
    tf_targets = c("NR3C1", "NR3C2"),
    legend_title = "Regulation of GR/MR-dependent genes",
    y_range = NULL
) {
  # Pomocnicza funkcja: zliczanie unikalnych genów dla TF
  count_genes_by_tf <- function(df, tf_name) {
    df %>%
      filter(
        Adjusted.P.value < 0.05,
        n_genes > 2,
        grepl(tf_name, Term, ignore.case = TRUE)
      ) %>%
      pull(Genes) %>%
      strsplit(";") %>%
      unlist() %>%
      unique() %>%
      length()
  }
  
  # Tworzenie podsumowania
  summary_df <- map_dfr(names(enrichr_results), function(tissue) {
    res <- enrichr_results[[tissue]]
    expand.grid(
      regulation = names(res),
      tf = tf_targets,
      stringsAsFactors = FALSE
    ) %>%
      rowwise() %>%
      mutate(
        tissue = tissue,
        n_genes = count_genes_by_tf(res[[regulation]], tf),
        x_label = tf
      )
  }) %>%
    ungroup()
  
  # Uporządkowanie
  summary_df$tissue <- factor(summary_df$tissue, levels = c("NeuralCells", "LungCells", "BloodCells"))
  summary_df$x_label <- factor(summary_df$x_label, levels = tf_targets)
  summary_df$regulation <- factor(summary_df$regulation, levels = c("Up", "Down"))
  
  # Sumy dla etykiet
  summary_sum <- summary_df %>%
    group_by(tissue, x_label) %>%
    summarise(
      max_genes = max(n_genes),
      total_genes = sum(n_genes),
      .groups = "drop"
    )
  
  # Wykres
  p <- ggplot(summary_df, aes(x = x_label, y = n_genes, fill = regulation)) +
    geom_col(position = position_dodge(width = 0.7), color = "black", width = 0.7) +
    geom_text(
      data = summary_sum,
      aes(x = x_label, y = max_genes + 6, label = total_genes),
      inherit.aes = FALSE,
      fontface = "bold",
      size = 4,
      color = "black"
    ) +
    facet_wrap(~ tissue, nrow = 1) +
    theme_classic(base_size = 14) +
    labs(
      x = paste("Aggregated", unique(summary_df$tissue), "terms related to", paste(tf_targets, collapse = " / ")),
      y = "Number of unique genes",
      fill = legend_title
    ) +
    scale_fill_manual(values = c("Up" = "firebrick", "Down" = "darkblue")) +
    theme(
      text = element_text(color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black", face = "bold"),
      strip.text = element_text(face = "bold", size = 14, color = "black"),
      legend.position = "bottom",
      legend.title = element_text(size = 12, face = "bold", color = "black"),
      legend.text = element_text(size = 11, color = "black"),
      axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
      axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 10)),
      axis.title.y = element_text(size = 13),
      panel.border = element_rect(color = "black", fill = NA, size = 0.8)
    )
  
  # Zakres osi Y (opcjonalny)
  if (!is.null(y_range)) {
    p <- p + coord_cartesian(ylim = y_range)
  }
  
  return(p)
}

# ============================================================
# 🧬 2. Uruchomienie Enrichr dla wybranych tkanek
# ============================================================
tissues <- c("NeuralCells", "LungCells", "BloodCells")

# Przykład dla bazy ChEA_2022 (bez prefiksów)
enrichr_chea_fullSignatures <- tissues %>%
  set_names() %>%
  map(~ get_enrichr_results_generic(
    tissue = .x,
    database = "ChEA_2022"
  ))

enrichr_chea_minusGlobal6TissueDerived <- tissues %>%
  set_names() %>%
  map(~ get_enrichr_results_generic(
    tissue = .x,
    database = "ChEA_2022",
    prefix_up = "minusGlobalUp6TissuesDerivedCells_",
    prefix_down = "minusGlobalDown6TissuesDerivedCells_"
  ))

enrichr_chea_minusGlobal5TissueDerived <- tissues %>%
  set_names() %>%
  map(~ get_enrichr_results_generic(
    tissue = .x,
    database = "ChEA_2022",
    prefix_up = "minusGlobalUp5TissuesDerivedCells_",
    prefix_down = "minusGlobalDown5TissuesDerivedCells_"
  ))

enrichr_chea_minusClusterBMC <- tissues %>%
  set_names() %>%
  map(~ get_enrichr_results_generic(
    tissue = .x,
    database = "ChEA_2022",
    prefix_up = "minusClustersKPO_",
    prefix_down = "minusClusterD_"
  ))

# ============================================================
# 📊 4. Przykład użycia – ChEA_2022, NR3C1/NR3C2
# ============================================================
p1 <- plot_enrichr_tf_results(
  enrichr_results = enrichr_chea_fullSignatures,
  tf_targets = c("NR3C1", "NR3C2"),
  legend_title = "Regulation of GR/MR-dependent genes",
  y_range = NULL
) 
p2 <- plot_enrichr_tf_results(
  enrichr_results = enrichr_chea_minusClusterBMC,
  tf_targets = c("NR3C1", "NR3C2"),
  legend_title = "Regulation of GR/MR-dependent genes",
  y_range = NULL
)

p3 <- plot_enrichr_tf_results(
  enrichr_results = enrichr_chea_minusGlobal6TissueDerived,
  tf_targets = c("NR3C1", "NR3C2"),
  legend_title = "Regulation of GR/MR-dependent genes",
  y_range = NULL
)

p4 <- plot_enrichr_tf_results(
  enrichr_results = enrichr_chea_minusGlobal5TissueDerived,
  tf_targets = c("NR3C1", "NR3C2"),
  legend_title = "Regulation of GR/MR-dependent genes",
  y_range = NULL
)

p1 + p2 + p3 + p4


# ============================================================
# 🏷️ Dodaj proste robocze tytuły do każdego wykresu
# ============================================================
p1 <- p1 + ggtitle("Full signatures")
p2 <- p2 + ggtitle("Minus cluster BMC")
p3 <- p3 + ggtitle("Minus global (6 tissues)")
p4 <- p4 + ggtitle("Minus global (5 tissues)")

# ============================================================
# 🧩 Połączenie wszystkich wykresów z jedną legendą na dole
# ============================================================
combined_plot <- ((p1 | p2) / (p3 | p4)) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "ChEA_2022 enrichment for NR3C1 / NR3C2") &
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(size = 12),
    plot.caption = element_text(size = 10),
    plot.title.position = "panel"
  )

# Wyświetlenie
combined_plot



# ============================================================
# 📂 Ścieżki i ustawienia
# ============================================================
output_dir <- "/home/mateusz/projects/ifpan-michkor-gr/results_v2/positive_validation_grSignatures_21.10.2025/chea_cellmaker"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

today <- "21.10.2025"

# ============================================================
# 💾 Zapis wyników Enrichr (RDS)
# ============================================================
saveRDS(enrichr_chea_fullSignatures,
        file.path(output_dir, paste0("tissueCellsGrSignatures_chea_fullSignatures_", today, ".rds")))
saveRDS(enrichr_chea_minusClusterBMC,
        file.path(output_dir, paste0("tissueCellsGrSignatures_chea_minusClusterBMC_", today, ".rds")))
saveRDS(enrichr_chea_minusGlobal6TissueDerived,
        file.path(output_dir, paste0("tissueCellsGrSignatures_chea_minusGlobal6TissueDerived_", today, ".rds")))
saveRDS(enrichr_chea_minusGlobal5TissueDerived,
        file.path(output_dir, paste0("tissueCellsGrSignatures_chea_minusGlobal5TissueDerived_", today, ".rds")))

# ============================================================
# 💾 Zapis wykresu (SVG)
# ============================================================
svg_filename <- file.path(output_dir, paste0("tissueCellsGrSignatures_chea_combinedPlot_", today, ".svg"))
svg(svg_filename, width = 14, height = 10)
print(combined_plot)
dev.off()

message("✅ Wyniki zapisane w katalogu: ", output_dir)




# ============================================================
# 📊 5. Wersja z ustaloną osią Y (y_range = c(0, 220))
# ============================================================

p1_yFreeze <- plot_enrichr_tf_results(
  enrichr_results = enrichr_chea_fullSignatures,
  tf_targets = c("NR3C1", "NR3C2"),
  legend_title = "Regulation of GR/MR-dependent genes",
  y_range = c(0, 220)
) + ggtitle("Full signatures")

p2_yFreeze <- plot_enrichr_tf_results(
  enrichr_results = enrichr_chea_minusClusterBMC,
  tf_targets = c("NR3C1", "NR3C2"),
  legend_title = "Regulation of GR/MR-dependent genes",
  y_range = c(0, 220)
) + ggtitle("Minus cluster BMC")

p3_yFreeze <- plot_enrichr_tf_results(
  enrichr_results = enrichr_chea_minusGlobal6TissueDerived,
  tf_targets = c("NR3C1", "NR3C2"),
  legend_title = "Regulation of GR/MR-dependent genes",
  y_range = c(0, 220)
) + ggtitle("Minus global (6 tissues)")

p4_yFreeze <- plot_enrichr_tf_results(
  enrichr_results = enrichr_chea_minusGlobal5TissueDerived,
  tf_targets = c("NR3C1", "NR3C2"),
  legend_title = "Regulation of GR/MR-dependent genes",
  y_range = c(0, 220)
) + ggtitle("Minus global (5 tissues)")

# ============================================================
# 🧩 Połączenie z jedną legendą
# ============================================================
combined_plot_yFreeze <- ((p1_yFreeze | p2_yFreeze) / (p3_yFreeze | p4_yFreeze)) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "ChEA_2022 enrichment for NR3C1 / NR3C2 (Y fixed 0–220)") &
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(size = 12),
    plot.caption = element_text(size = 10),
    plot.title.position = "panel"
  )

# ============================================================
# 💾 Zapis SVG (yFreeze)
# ============================================================
output_dir <- "/home/mateusz/projects/ifpan-michkor-gr/results_v2/positive_validation_grSignatures_21.10.2025/chea_cellmaker"
today <- "21.10.2025"

svg_filename_yFreeze <- file.path(output_dir, paste0("tissueCellsGrSignatures_chea_combinedPlot_yFreeze_", today, ".svg"))
svg(svg_filename_yFreeze, width = 14, height = 10)
print(combined_plot_yFreeze)
dev.off()

message("✅ Wykres z ustaloną osią Y zapisany jako: ", svg_filename_yFreeze)
