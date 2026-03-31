
# ============================================================
# 🧩 Definicja tkanek i nazw
# ============================================================
tissues <- c("BloodCells", "NeuralCells", "LungCells")
term_categories <- c("Brain", "Lung", "Blood")  # szukane słowa w kolumnie Term

# ============================================================
# 🧠 Funkcja pomocnicza do tworzenia listy wzbogaceń
# ============================================================
get_enrichr_results <- function(tissue) {
  up_list <- flat_allGrSignatures_17.10.2025[[paste0("minusGlobalUp5TissuesDerivedCells_", tissue, "Up")]]
  down_list <- flat_allGrSignatures_17.10.2025[[paste0("minusGlobalDown5TissuesDerivedCells_", tissue, "Down")]]
  
  list(
    Up = run_enrichr(gene_list = up_list, database = "CellMarker_2024"),
    Down = run_enrichr(gene_list = down_list, database = "CellMarker_2024"),
    Both = run_enrichr(gene_list = c(up_list, down_list), database = "CellMarker_2024")
  )
}

# ============================================================
# 📊 Wykonanie analiz wzbogacenia dla wszystkich tkanek
# ============================================================
enrichr_all <- tissues %>%
  set_names() %>%
  map(get_enrichr_results)

# ============================================================
# 🧮 Funkcja do zliczania liczby unikalnych genów dla kategorii
# ============================================================
count_genes_by_category <- function(df, keyword) {
  df %>%
    filter(
      Adjusted.P.value < 0.05,
      n_genes > 2,
      grepl(keyword, Term, ignore.case = TRUE)
    ) %>%
    pull(Genes) %>%
    strsplit(";") %>%
    unlist() %>%
    unique() %>%
    length()
}

# ============================================================
# 📈 Tworzenie podsumowania
# ============================================================
summary_df <- map_dfr(names(enrichr_all), function(tissue) {
  res <- enrichr_all[[tissue]]
  expand.grid(
    regulation = names(res),
    category = term_categories,
    stringsAsFactors = FALSE
  ) %>%
    rowwise() %>%
    mutate(
      tissue = tissue,
      n_genes = count_genes_by_category(res[[regulation]], category)
    )
}) %>%
  ungroup()

# ============================================================
# 🎨 Wykres słupkowy
# ============================================================
ggplot(summary_df, aes(x = tissue, y = n_genes, fill = category)) +
  geom_col(position = "dodge", color = "black") +
  facet_wrap(~ regulation, ncol = 1) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Unique Enrichr CellMarker_2024 genes by tissue and category",
    subtitle = "Filtered for Adjusted.P.value < 0.05 and n_genes > 2",
    x = "Tissue (source of GR-dependent signature)",
    y = "Number of unique genes",
    fill = "Cell-type keyword"
  ) +
  scale_fill_manual(values = c("Brain" = "#d73027", "Lung" = "#4575b4", "Blood" = "#91bfdb")) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )
