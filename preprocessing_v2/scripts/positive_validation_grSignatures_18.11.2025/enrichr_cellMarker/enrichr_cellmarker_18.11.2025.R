# ============================================================
# 📦 Wczytanie pakietów
# ============================================================
library(dplyr)
library(purrr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(openxlsx)

tissues <- c("NeuralCells", "LungCells", "BloodCells")

# ============================================================
# 🧠 Pobieranie list genów zgodnie z REALNYMI nazwami w flatten
# ============================================================
get_genes <- function(prefix, tissue, type) {
  flat_allGrSignatures_18.11.2025[[paste0(prefix, tissue, type)]]
}

# ============================================================
# 🧪 Przygotowanie wyników do Enrichr
# ============================================================
run_variant <- function(prefix) {
  map(tissues, function(tissue) {
    list(
      Up   = get_genes(prefix, tissue, "Up"),
      Down = get_genes(prefix, tissue, "Down")
    )
  }) %>% set_names(tissues)
}

# --- ZESTAWY LIST GENÓW --- #
geneLists_raw   <- run_variant("")
geneLists_minus4 <- run_variant("minusGlobalUpDown4TissuesDerivedCells_")
geneLists_minus5 <- run_variant("minusGlobalUpDown5TissuesDerivedCells_")
geneLists_minus6 <- run_variant("minusGlobalUpDown6TissuesDerivedCells_")

# ============================================================
# ENRICHR
# ============================================================

run_cellmarker <- function(list_set) {
  map(list_set, ~ list(
    Up   = run_enrichr(.x$Up,   database = "CellMarker_2024"),
    Down = run_enrichr(.x$Down, database = "CellMarker_2024")
  ))
}

enrichr_raw   <- run_cellmarker(geneLists_raw)
enrichr_minus4 <- run_cellmarker(geneLists_minus4)
enrichr_minus5 <- run_cellmarker(geneLists_minus5)
enrichr_minus6 <- run_cellmarker(geneLists_minus6)

# ============================================================
# 📊 WYKRES CELL MARKER
# ============================================================
plot_enrichr_cellmarker_results <- function(
    enrichr_results,
    term_categories = c("Brain", "Lung", "Blood"),
    legend_title = "Regulation of GR-dependent genes",
    y_range = NULL,
    min_overlap_threshold = 2,
    fdr_threshold = 0.05
) {
  # pomocnicza funkcja do zliczania unikalnych genów na kategorię
  count_genes_by_category <- function(df, keyword) {
    df %>%
      filter(
        Adjusted.P.value < fdr_threshold,
        n_genes > min_overlap_threshold,
        grepl(keyword, Term, ignore.case = TRUE)
      ) %>%
      pull(Genes) %>%
      strsplit(";") %>%
      unlist() %>%
      trimws() %>%
      unique() %>%
      length()
  }
  
  # Tworzenie podsumowania
  summary_df <- purrr::map_dfr(names(enrichr_results), function(tissue) {
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
  
  # Sumy do etykiet
  summary_sum <- summary_df %>%
    group_by(tissue, x_label) %>%
    summarise(
      max_genes = max(n_genes),
      total_genes = sum(n_genes),
      .groups = "drop"
    )
  
  # wykres
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


p_raw   <- plot_enrichr_cellmarker_results(enrichr_raw, fdr_threshold = 0.05, min_overlap_threshold = 3)   + ggtitle("RAW GR-dependent signatures")
p_4     <- plot_enrichr_cellmarker_results(enrichr_minus4, fdr_threshold = 0.05, min_overlap_threshold = 3) + ggtitle("MinusGlobalUpDown (4 tissues)")
p_5     <- plot_enrichr_cellmarker_results(enrichr_minus5, fdr_threshold = 0.1, min_overlap_threshold = 3) + ggtitle("MinusGlobalUpDown (5 tissues)")
p_6     <- plot_enrichr_cellmarker_results(enrichr_minus6, fdr_threshold = 0.05, min_overlap_threshold = 3) + ggtitle("MinusGlobalUpDown (6 tissues)")

combined_plot <- wrap_plots(p_raw, p_4, p_5, p_6, ncol = 2)

# ============================================================
# ZAPIS SVG i RDS
# ============================================================

output_dir <- "/home/mateusz/projects/ifpan-michkor-gr/results_v2/positive_validation_grSignatures_18.11.2025/enrichr_cellmarker"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

svg_file <- file.path(output_dir, "cellmarker_minusGlobalUpDown_RAW_4_5_6_18.11.2025.svg")
svg(svg_file, width = 14, height = 10)
print(combined_plot)
dev.off()


# ============================================================
# WERSJA DLA FDR 0.01
# ============================================================

p_raw_fdr01 <- plot_enrichr_cellmarker_results(
  enrichr_raw,
  fdr_threshold = 0.01,
  min_overlap_threshold = 3
) + ggtitle("RAW GR-dependent signatures (FDR < 0.01)")

p_4_fdr01 <- plot_enrichr_cellmarker_results(
  enrichr_minus4,
  fdr_threshold = 0.01,
  min_overlap_threshold = 3
) + ggtitle("MinusGlobalUpDown (4 tissues) (FDR < 0.01)")

p_5_fdr01 <- plot_enrichr_cellmarker_results(
  enrichr_minus5,
  fdr_threshold = 0.01,
  min_overlap_threshold = 3
) + ggtitle("MinusGlobalUpDown (5 tissues) (FDR < 0.01)")

p_6_fdr01 <- plot_enrichr_cellmarker_results(
  enrichr_minus6,
  fdr_threshold = 0.01,
  min_overlap_threshold = 3
) + ggtitle("MinusGlobalUpDown (6 tissues) (FDR < 0.01)")

combined_plot_fdr01 <- wrap_plots(
  p_raw_fdr01, p_4_fdr01, p_5_fdr01, p_6_fdr01,
  ncol = 2
)

# ============================================================
# ZAPIS SVG
# ============================================================

svg_file_fdr01 <- file.path(
  output_dir,
  "cellmarker_minusGlobalUpDown_RAW_4_5_6_FDR01_18.11.2025.svg"
)

svg(svg_file_fdr01, width = 14, height = 10)
print(combined_plot_fdr01)
dev.off()

saveRDS(list(
  RAW   = enrichr_raw,
  minus4 = enrichr_minus4,
  minus5 = enrichr_minus5,
  minus6 = enrichr_minus6
), file.path(output_dir, "cellmarker_allVariants_18.11.2025.rds"))


# ===============================================
# 4 pliki XLSX
# ===============================================
format_enrichr_df <- function(df) {
  df %>%
    # wyciąganie liczby przed slashem
    mutate(
      n_genes_overlap = as.numeric(sub("/.*", "", Overlap))
    ) %>%
    # filtr
    filter(
      n_genes_overlap >= 2,
      Adjusted.P.value < 0.05
    ) %>%
    # usuń stare pola
    select(-Old.P.value, -Old.Adjusted.P.value) %>%
    # formatowanie wartości
    mutate(
      Odds.Ratio       = round(Odds.Ratio, 3),
      Combined.Score   = round(Combined.Score, 3),
      P.value          = formatC(P.value, format = "e", digits = 3),
      Adjusted.P.value = formatC(Adjusted.P.value, format = "e", digits = 3)
    )
}

output_dir <- "/home/mateusz/projects/ifpan-michkor-gr/results_v2/positive_validation_grSignatures_18.11.2025/enrichr_cellmarker"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

result_list <- list(
  RAW    = enrichr_raw,
  minus4 = enrichr_minus4,
  minus5 = enrichr_minus5,
  minus6 = enrichr_minus6
)

walk(names(result_list), function(variant) {
  
  wb <- createWorkbook()
  obj <- result_list[[variant]]
  
  iwalk(obj, function(tissue_obj, tissue_name) {
    iwalk(tissue_obj, function(df, regulation_name) {
      
      df <- format_enrichr_df(df)  # FORMAT + FILTER
      
      sheet_name <- paste(tissue_name, regulation_name, sep = "_")
      sheet_name <- str_sub(sheet_name, 1, 31)
      
      addWorksheet(wb, sheet_name)
      writeData(wb, sheet_name, df, withFilter = TRUE)
      freezePane(wb, sheet = sheet_name, firstActiveRow = 2)
      setColWidths(wb, sheet = sheet_name, cols = 1:ncol(df), widths = "auto")
    })
  })
  
  saveWorkbook(
    wb,
    file = file.path(output_dir, paste0("cellmarker_", variant, "_18.11.2025.xlsx")),
    overwrite = TRUE
  )
})

