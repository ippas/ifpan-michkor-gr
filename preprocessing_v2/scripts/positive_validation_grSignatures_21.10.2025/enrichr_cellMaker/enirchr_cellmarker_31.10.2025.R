# ============================================================
# 📦 Wczytanie pakietów
# ============================================================
library(dplyr)
library(purrr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(openxlsx)

# ============================================================
# 🧠 1. Uniwersalna funkcja do analizy Enrichr (CellMarker)
# ============================================================
get_enrichr_results_generic <- function(
    tissue,
    database = "CellMarker_2024",
    prefix_up = "",
    prefix_down = ""
) {
  up_name <- paste0(prefix_up, tissue, "Up")
  down_name <- paste0(prefix_down, tissue, "Down")
  
  up_list <- flat_allGrSignatures_31.10.2025[[up_name]]
  down_list <- flat_allGrSignatures_31.10.2025[[down_name]]
  
  if (is.null(up_list)) warning(paste("⚠️ Nie znaleziono listy:", up_name))
  if (is.null(down_list)) warning(paste("⚠️ Nie znaleziono listy:", down_name))
  
  list(
    Up = run_enrichr(gene_list = up_list, database = database),
    Down = run_enrichr(gene_list = down_list, database = database)
  )
}

# ============================================================
# 🧩 2. Funkcja do tworzenia wykresów CellMarker (Brain/Lung/Blood)
# ============================================================
plot_enrichr_cellmarker_results <- function(
    enrichr_results,
    term_categories = c("Brain", "Lung", "Blood"),
    legend_title = "Regulation of GR-dependent genes",
    y_range = NULL
) {
  # Pomocnicza funkcja do zliczania unikalnych genów na kategorię
  count_genes_by_category <- function(df, keyword) {
    df %>%
      filter(Adjusted.P.value < 0.05, n_genes > 2, grepl(keyword, Term, ignore.case = TRUE)) %>%
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
      category = term_categories,
      stringsAsFactors = FALSE
    ) %>%
      rowwise() %>%
      mutate(
        tissue = tissue,
        n_genes = count_genes_by_category(res[[regulation]], category),
        x_label = category
      )
  }) %>% ungroup()
  
  summary_df$tissue <- factor(summary_df$tissue, levels = c("NeuralCells", "LungCells", "BloodCells"))
  summary_df$x_label <- factor(summary_df$x_label, levels = term_categories)
  summary_df$regulation <- factor(summary_df$regulation, levels = c("Up", "Down"))
  
  # Sumy dla etykiet
  summary_sum <- summary_df %>%
    group_by(tissue, x_label) %>%
    summarise(max_genes = max(n_genes), total_genes = sum(n_genes), .groups = "drop")
  
  # Wykres
  p <- ggplot(summary_df, aes(x = x_label, y = n_genes, fill = regulation)) +
    geom_col(position = position_dodge(width = 0.7), color = "black", width = 0.7) +
    geom_text(
      data = summary_sum,
      aes(x = x_label, y = max_genes + 2, label = total_genes),
      inherit.aes = FALSE,
      fontface = "bold",
      size = 4,
      color = "black"
    ) +
    facet_wrap(~tissue, nrow = 1) +
    theme_classic(base_size = 14) +
    labs(
      x = "Aggregated results for tissue-related CellMarker terms",
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
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
      axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 10)),
      axis.title.y = element_text(size = 13),
      panel.border = element_rect(color = "black", fill = NA, size = 0.8)
    )
  
  if (!is.null(y_range)) {
    p <- p + coord_cartesian(ylim = y_range)
  }
  
  return(p)
}

# ============================================================
# 🧬 3. Uruchomienie Enrichr dla CellMarker_2024 (Minus global UP+DOWN)
# ============================================================

tissues <- c("NeuralCells", "LungCells", "BloodCells")

enrichr_cellmarker_minusGlobalUpDown5TissueDerived <- tissues %>%
  set_names() %>%
  map(~ get_enrichr_results_generic(
    tissue = .x,
    database = "CellMarker_2024",
    prefix_up = "minusGlobalUpDown5TissuesDerivedCells_",
    prefix_down = "minusGlobalUpDown5TissuesDerivedCells_"
  ))

# ============================================================
# 📊 4. Tworzenie wykresu
# ============================================================

p_updown <- plot_enrichr_cellmarker_results(enrichr_cellmarker_minusGlobalUpDown5TissueDerived) +
  ggtitle("Minus global (UP+DOWN, 5 tissues)")

# ============================================================
# 📂 Zapisywanie wyników
# ============================================================

output_dir <- "/home/mateusz/projects/ifpan-michkor-gr/results_v2/positive_validation_grSignatures_31.10.2025/enrichr_cellmarker_updown"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
today <- "31.10.2025"

# 💾 Zapis RDS
saveRDS(enrichr_cellmarker_minusGlobalUpDown5TissueDerived,
        file.path(output_dir, paste0("tissueCellsGrSignatures_cellmarker_minusGlobalUpDown5TissueDerived_", today, ".rds")))

# 💾 Zapis SVG
svg_filename <- file.path(output_dir, paste0("tissueCellsGrSignatures_cellmarker_minusGlobalUpDown5TissueDerived_plot_", today, ".svg"))
svg(svg_filename, width = 8, height = 6)
print(p_updown)
dev.off()
message("✅ Wykres zapisany: ", svg_filename)

# ============================================================
# 📘 5. Eksport do XLSX (wyniki z enrichr)
# ============================================================

filtered_list <- list(
  "NeuralCellsUp"   = enrichr_cellmarker_minusGlobalUpDown5TissueDerived$NeuralCells$Up,
  "NeuralCellsDown" = enrichr_cellmarker_minusGlobalUpDown5TissueDerived$NeuralCells$Down,
  "BloodCellsUp"    = enrichr_cellmarker_minusGlobalUpDown5TissueDerived$BloodCells$Up,
  "BloodCellsDown"  = enrichr_cellmarker_minusGlobalUpDown5TissueDerived$BloodCells$Down,
  "LungCellsUp"     = enrichr_cellmarker_minusGlobalUpDown5TissueDerived$LungCells$Up,
  "LungCellsDown"   = enrichr_cellmarker_minusGlobalUpDown5TissueDerived$LungCells$Down
) %>%
  lapply(function(x) {
    x %>%
      filter(n_genes > 2, Adjusted.P.value < 0.05) %>%
      select(-c(Old.P.value, Old.Adjusted.P.value))
  })

output_file <- file.path(
  output_dir,
  "cellMarker2024_grSignaturesMinusGlobalUpDown5TissuesDerived_fdr0.05overlap3-31.10.2025.xlsx"
)

wb <- createWorkbook()
sheet_order <- names(filtered_list)
style_3dec <- createStyle(numFmt = "0.000")

for (sheet in sheet_order) {
  df <- filtered_list[[sheet]]
  addWorksheet(wb, sheet)
  writeData(wb, sheet, df, withFilter = TRUE)
  freezePane(wb, sheet, firstRow = TRUE)
  num_cols <- intersect(c("P.value", "Adjusted.P.value", "Odds.Ratio", "Combined.Score"), colnames(df))
  for (col_name in num_cols) {
    col_index <- which(colnames(df) == col_name)
    addStyle(wb, sheet, style = style_3dec, rows = 2:(nrow(df) + 1), cols = col_index, gridExpand = TRUE)
  }
}

saveWorkbook(wb, output_file, overwrite = TRUE)
message("✅ Zapisano plik XLSX: ", output_file)
