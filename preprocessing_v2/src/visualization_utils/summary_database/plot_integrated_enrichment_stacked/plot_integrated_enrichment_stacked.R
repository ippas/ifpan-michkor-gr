plot_integrated_enrichment_stacked <- function(
    gene_summary_subset,
    subdata = "summary_up",
    enrichr_source = "enrichr_3pub",
    fdr_threshold = 0.05,
    y_limits = c(0, 60),
    palette_sources = c(
      "ChEA_2022" = "#4b0b0b",
      "LINCS_L1000" = "#9c5757",
      "CellMarker_2024" = "#f0c1c1"
    ),
    # --- parametry estetyczne ---
    axis_text_size = 12,   # pt
    axis_title_size = 13,  # pt
    label_text_size = 10,  # pt (konwertowany)
    legend_text_size = 10, # pt
    y_axis_label = expression(-log[10]("FDR")),
    x_text_angle = 15      # stopnie nachylenia etykiet osi X
) {
  library(dplyr)
  library(ggplot2)
  library(stringr)
  
  # --- funkcja pomocnicza: konwersja punktów (pt) → jednostki ggplot2 ---
  pt_to_gg_size <- function(pt) pt / 2.845
  
  # --- 1️⃣ Pobranie danych enrichr ---
  enrich_data <- gene_summary_subset[[subdata]]$enrichr_summary_source[[enrichr_source]]
  if (is.null(enrich_data)) {
    stop(paste0("❌ Nie znaleziono źródła: ", enrichr_source))
  }
  
  # --- 2️⃣ Kierunek (Up / Down) ---
  direction <- ifelse(grepl("down", subdata, ignore.case = TRUE), "Down", "Up")
  
  # --- 3️⃣ ChEA ---
  df_chea <- tryCatch({
    enrich_data$ChEA_2022 %>%
      filter(Adjusted.P.value < fdr_threshold) %>%
      filter(grepl("NR3C1|NR3C2", Term, ignore.case = TRUE)) %>%
      mutate(
        category = case_when(
          str_detect(Term, "NR3C1") ~ "NR3C1",
          str_detect(Term, "NR3C2") ~ "NR3C2"
        ),
        source = "ChEA_2022",
        log10FDR = -log10(Adjusted.P.value)
      )
  }, error = function(e) data.frame())
  
  # --- 4️⃣ LINCS ---
  lincs_terms <- paste(
    c("Methylprednisolone", "Dexamethasone", "Budesonide", "Hydrocortisone"),
    direction
  )
  
  df_lincs <- tryCatch({
    enrich_data$LINCS_L1000_Chem_Pert_Consensus_Sigs %>%
      filter(Adjusted.P.value < fdr_threshold) %>%
      filter(Term %in% lincs_terms) %>%
      mutate(
        category = Term,
        source = "LINCS_L1000",
        log10FDR = -log10(Adjusted.P.value)
      )
  }, error = function(e) data.frame())
  
  # --- 5️⃣ CellMarker ---
  df_cell <- tryCatch({
    enrich_data$CellMarker_2024 %>%
      filter(Adjusted.P.value < fdr_threshold) %>%
      filter(grepl("Astrocyte|Neuron|Microglial Cell|Oligodendrocyte", Term, ignore.case = TRUE)) %>%
      filter(grepl("Brain", Term, ignore.case = TRUE)) %>%
      mutate(
        category = case_when(
          str_detect(Term, "Astrocyte") ~ "Astrocyte",
          str_detect(Term, "Neuron") ~ "Neuron",
          str_detect(Term, "Microglial Cell") ~ "Microglia",
          str_detect(Term, "Oligodendrocyte Precursor Cell") ~ "OPC",
          str_detect(Term, "Oligodendrocyte") ~ "Oligodendrocyte"
        ),
        source = "CellMarker_2024",
        log10FDR = -log10(Adjusted.P.value)
      )
  }, error = function(e) data.frame())
  
  # --- 6️⃣ Połączenie danych ---
  df_plot_raw <- bind_rows(df_chea, df_lincs, df_cell)
  if (nrow(df_plot_raw) == 0) stop("❌ Brak danych do wyświetlenia po filtracji FDR.")
  
  # --- 7️⃣ Kategorie ---
  cat_levels <- c(
    "NR3C1", "NR3C2",
    paste("Methylprednisolone", direction),
    paste("Dexamethasone", direction),
    paste("Budesonide", direction),
    paste("Hydrocortisone", direction),
    "Astrocyte", "Neuron", "Microglia", "Oligodendrocyte", "OPC"
  )
  
  # --- 8️⃣ Agregacja ---
  df_plot <- df_plot_raw %>%
    group_by(source, category) %>%
    summarise(
      log10FDR = max(log10FDR, na.rm = TRUE),
      n_terms = n(),
      n_genes = length(unique(unlist(str_split(Genes, ";")))),
      .groups = "drop"
    )
  
  # --- 9️⃣ Uzupełnienie braków ---
  all_combinations <- expand.grid(
    source = names(palette_sources),
    category = cat_levels,
    stringsAsFactors = FALSE
  )
  
  df_plot_full <- all_combinations %>%
    left_join(df_plot, by = c("source", "category")) %>%
    mutate(
      log10FDR = ifelse(is.na(log10FDR), 0, log10FDR),
      n_genes = ifelse(is.na(n_genes), 0, n_genes),
      n_terms = ifelse(is.na(n_terms), 0, n_terms)
    )
  
  # --- 🔟 Etykiety ---
  df_labels <- df_plot_full %>%
    group_by(category) %>%
    summarise(
      total_genes = sum(n_genes),
      total_terms = sum(n_terms),
      total_height = sum(log10FDR)
    ) %>%
    mutate(
      label = ifelse(
        total_terms > 1,
        paste0(total_genes, " genes\n(", total_terms, " terms)"),
        paste0(total_genes, " genes")
      )
    )
  
  df_plot_full$category <- factor(df_plot_full$category, levels = cat_levels)
  
  # --- 11️⃣ Wykres ---
  p <- ggplot(df_plot_full, aes(x = category, y = log10FDR, fill = source)) +
    geom_col(color = "black", width = 0.8) +
    geom_text(
      data = df_labels,
      aes(x = category, y = total_height, label = label),
      vjust = -0.5,
      size = pt_to_gg_size(label_text_size),
      color = "black",                    # ✅ czarny tekst nad słupkami
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = palette_sources, drop = FALSE) +
    labs(
      x = NULL,
      y = y_axis_label
    ) +
    theme_classic(base_size = axis_text_size) +
    theme(
      axis.text.x = element_text(
        face = "bold",
        angle = x_text_angle,
        hjust = ifelse(x_text_angle == 0, 0.5, 1),
        size = axis_text_size,
        color = "black"                    # ✅ czarna oś X
      ),
      axis.text.y = element_text(size = axis_text_size, color = "black"),  # ✅ czarna oś Y
      axis.title.y = element_text(size = axis_title_size, face = "bold", color = "black"), # ✅ czarny tytuł osi
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = legend_text_size, color = "black"), # ✅ czarna legenda
      plot.title = element_text(color = "black"),
      plot.subtitle = element_text(color = "black"),
      plot.caption = element_text(color = "black")
    ) +
    coord_cartesian(clip = "off", ylim = y_limits)
  
  return(p)
}
