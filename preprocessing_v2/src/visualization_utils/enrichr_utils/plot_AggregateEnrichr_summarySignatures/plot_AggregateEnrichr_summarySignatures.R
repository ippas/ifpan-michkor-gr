plot_AggregateEnrichr_summarySignatures <- function(
    data,
    subdata = "summary_down",
    gene_set_name = "genes_3pub",      # 👈 nowy argument
    enrichr_set_name = "enrichr_3pub", # 👈 nowy argument
    plot_title = NULL,
    x_label = NULL,
    y_label = NULL,
    palette_sources = c(
      "ChEA_2022" = "#4b0b0b",
      "LINCS_L1000" = "#9c5757",
      "CellMarker_2024" = "#f0c1c1"
    ),
    axis_text_size = 14,
    axis_title_size = 16,
    legend_text_size = 14,
    title_size = 18,
    signature_text_size = 5,
    number_association = TRUE,
    y_limits = NULL
) {
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(purrr)
  
  signature_text_size <- signature_text_size / 2.835
  
  # --- 1️⃣ liczba genów w sygnaturze (dynamicznie)
  if (!gene_set_name %in% names(data[[subdata]]$enrichr_summary_source)) {
    stop(paste0("Gene set '", gene_set_name, "' not found in enrichr_summary_source"))
  }
  n_signature_genes <- length(data[[subdata]]$enrichr_summary_source[[gene_set_name]])
  
  # --- 2️⃣ sprawdź, czy enrichr_set_name istnieje
  if (!enrichr_set_name %in% names(data[[subdata]]$enrichr_summary_source)) {
    stop(paste0("Enrichr set '", enrichr_set_name, "' not found in enrichr_summary_source"))
  }
  enrichr_data <- data[[subdata]]$enrichr_summary_source[[enrichr_set_name]]
  
  # --- 3️⃣ ChEA_2022 ---
  df_chee <- enrichr_data$ChEA_2022 %>%
    filter(Adjusted.P.value < 0.05, n_genes >= 3) %>%
    filter(grepl("NR3C1|NR3C2", Term, TRUE)) %>%
    mutate(
      aggregate_term = case_when(
        grepl("NR3C1", Term, TRUE) ~ "NR3C1",
        grepl("NR3C2", Term, TRUE) ~ "NR3C2",
        TRUE ~ "other"
      )
    ) %>%
    group_by(aggregate_term) %>%
    nest() %>%
    mutate(
      aggregate_genes    = map(data, ~ .x$Genes |> unlist() |> strsplit(";") |> unlist() |> unique()),
      n_aggreate_genes   = map_int(aggregate_genes, length),
      n_traits           = map_int(data, nrow),
      source             = "ChEA_2022"
    )
  
  # --- 4️⃣ LINCS_L1000 ---
  df_lincs <- enrichr_data$LINCS_L1000_Chem_Pert_Consensus_Sigs %>%
    filter(Adjusted.P.value < 0.05, n_genes >= 3) %>%
    filter(grepl("Betamethasone|Dexamethasone|Fluocortolone|Methylprednisolone|Paramethasone|Prednisolone|Prednisone|Triamcinolone|Hydrocortisone|Cortisone|Prednylidene|Rimexolone|Deflazacort|Cloprednol|Meprednisone|Cortivazol|Vamorolone", Term, TRUE)) %>%
    filter(grepl("Dexamethasone Up|Hydrocortisone Up|Methylprednisolone Up|Prednisone Up|Betamethasone Up|Budesonide Up|Dexamethasone Down|Hydrocortisone Down|Methylprednisolone Down|Prednisone Down|Betamethasone Up|Budesonide Down", Term, TRUE)) %>%
    mutate(
      aggregate_term = case_when(
        grepl("up", Term, TRUE) ~ "GC ↑",
        grepl("down", Term, TRUE) ~ "GC ↓",
        TRUE ~ "other"
      )
    ) %>%
    group_by(aggregate_term) %>%
    nest() %>%
    mutate(
      aggregate_genes    = map(data, ~ .x$Genes |> unlist() |> strsplit(";") |> unlist() |> unique()),
      n_aggreate_genes   = map_int(aggregate_genes, length),
      n_traits           = map_int(data, nrow),
      source             = "LINCS_L1000"
    )
  
  # --- 5️⃣ CellMarker_2024 ---
  df_cell <- enrichr_data$CellMarker_2024 %>%
    filter(Adjusted.P.value < 0.05, n_genes >= 3) %>%
    filter(grepl("lung|blood|brain", Term, TRUE)) %>%
    mutate(
      aggregate_term = case_when(
        grepl("brain", Term, TRUE) ~ "Brain cells",
        grepl("blood", Term, TRUE) ~ "Blood cells",
        grepl("lung",  Term, TRUE) ~ "Lung cells",
        TRUE ~ "other"
      )
    ) %>%
    group_by(aggregate_term) %>%
    nest() %>%
    mutate(
      aggregate_genes    = map(data, ~ .x$Genes |> unlist() |> strsplit(";") |> unlist() |> unique()),
      n_aggreate_genes   = map_int(aggregate_genes, length),
      n_traits           = map_int(data, nrow),
      source             = "CellMarker_2024"
    )
  
  # --- 6️⃣ połączenie
  combined_df <- bind_rows(df_chee, df_lincs, df_cell)
  
  # --- 7️⃣ uzupełnij brakujące terminy o 0
  term_levels <- c("NR3C1", "NR3C2", "GC ↑", "GC ↓", "Brain cells", "Blood cells", "Lung cells")
  missing_terms <- setdiff(term_levels, combined_df$aggregate_term)
  if (length(missing_terms) > 0) {
    missing_df <- tibble(
      aggregate_term = missing_terms,
      n_aggreate_genes = 0,
      n_traits = 0,
      source = c("ChEA_2022", "LINCS_L1000", "CellMarker_2024")[match(missing_terms, c("NR3C1", "NR3C2", "GC ↑", "GC ↓", "Brain cells", "Blood cells", "Lung cells"))]
    )
    combined_df <- bind_rows(combined_df, missing_df)
  }
  
  combined_df <- combined_df %>%
    mutate(aggregate_term = factor(aggregate_term, levels = term_levels))
  
  # --- 8️⃣ etykiety
  label_map <- combined_df %>%
    select(aggregate_term, n_traits) %>%
    distinct() %>%
    mutate(
      label = if (number_association)
        paste0(as.character(aggregate_term), " (n = ", n_traits, ")")
      else
        as.character(aggregate_term)
    ) %>%
    { setNames(.$label, .$aggregate_term) }
  
  # --- 9️⃣ wykres
  p <- ggplot(combined_df, aes(x = aggregate_term, y = n_aggreate_genes, fill = source)) +
    geom_bar(stat = "identity", width = 0.6, color = "black") +
    scale_x_discrete(labels = label_map) +
    scale_fill_manual(values = palette_sources) +
    theme_classic(base_size = axis_text_size) +
    theme(
      text = element_text(color = "black"),
      axis.text = element_text(color = "black", size = axis_text_size),
      axis.title = element_text(color = "black", size = axis_title_size),
      plot.title = element_text(color = "black", size = title_size, face = "bold"),
      legend.text = element_text(color = "black", size = legend_text_size),
      legend.title = element_blank(),
      legend.position = "bottom",
      axis.text.x = element_text(angle = 25, hjust = 1)
    )
  
  if (!is.null(y_limits)) p <- p + ylim(y_limits)
  
  # --- 🔟 liczba genów
  p <- p + annotate(
    "text",
    x = Inf, y = Inf,
    label = paste0("Number of genes in signature: ", n_signature_genes),
    hjust = 1.05, vjust = 1.5,
    size = signature_text_size,
    color = "black",
    fontface = "plain"
  )
  
  if (!is.null(plot_title)) p <- p + ggtitle(plot_title)
  if (!is.null(x_label))    p <- p + xlab(x_label)
  if (!is.null(y_label))    p <- p + ylab(y_label)
  
  return(p)
}

