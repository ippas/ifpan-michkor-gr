summarize_gene_dataset <- function(df, dataset_name = "dataset") {
  message("▶ Start summarize_gene_dataset for: ", dataset_name)
  
  df %>% ungroup() -> df
  stopifnot(is.data.frame(df))
  
  required_cols <- c("hgnc_symbol", "label", "source", "regulation")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("❌ Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # ------------------------------
  # Core summarization function
  # ------------------------------
  summarize_core <- function(data) {
    if (nrow(data) == 0) {
      return(list(
        n_records = 0, n_genes = 0, n_publications = 0, n_geneLists = 0,
        freq_genes_per_list = tibble(freq = integer(), n_genes = integer(), data = list()),
        freq_genes_per_paper = tibble(freq = integer(), n_genes = integer(), data = list())
      ))
    }
    
    result <- list(
      n_records      = nrow(data),
      n_genes        = n_distinct(data$hgnc_symbol),
      n_publications = n_distinct(data$source),
      n_geneLists    = n_distinct(data$label)
    )
    
    # ---- Liczba wystąpień genu w różnych listach ----
    freq_genes_per_list <- data %>%
      distinct(label, hgnc_symbol) %>%
      group_by(hgnc_symbol) %>%
      summarise(freq = n(), .groups = "drop") %>%
      group_by(freq) %>%
      summarise(
        n_genes = n(),
        data = list(hgnc_symbol),
        .groups = "drop"
      ) %>%
      arrange(freq)
    
    # ---- Liczba wystąpień genu w różnych publikacjach ----
    freq_genes_per_paper <- data %>%
      distinct(source, hgnc_symbol) %>%
      group_by(hgnc_symbol) %>%
      summarise(freq = n(), .groups = "drop") %>%
      group_by(freq) %>%
      summarise(
        n_genes = n(),
        data = list(hgnc_symbol),
        .groups = "drop"
      ) %>%
      arrange(freq)
    
    result$freq_genes_per_list  <- freq_genes_per_list
    result$freq_genes_per_paper <- freq_genes_per_paper
    return(result)
  }
  
  # ------------------------------
  # Regulation summary
  # ------------------------------
  regulation_summary <- df %>%
    filter(!is.na(regulation) & regulation %in% c("up", "down")) %>%
    distinct(hgnc_symbol, regulation) %>%
    group_by(hgnc_symbol) %>%
    summarise(
      has_up = any(regulation == "up"),
      has_down = any(regulation == "down"),
      regulation_type = case_when(
        has_up & has_down ~ "both",
        has_up & !has_down ~ "only_up",
        !has_up & has_down ~ "only_down",
        TRUE ~ NA_character_
      ),
      .groups = "drop"
    ) %>%
    group_by(regulation_type) %>%
    summarise(
      n_genes = n(),
      genes = list(hgnc_symbol),
      .groups = "drop"
    )
  
  # ------------------------------
  # Global summaries
  # ------------------------------
  summary_all  <- summarize_core(df)
  summary_up   <- summarize_core(df %>% filter(regulation == "up"))
  summary_down <- summarize_core(df %>% filter(regulation == "down"))
  
  # ------------------------------
  # Summary of gene counts per list (min, q1, median, mean, q3, max, sd)
  # ------------------------------
  n_genes_vec <- df %>%
    group_by(label) %>%
    summarise(n_genes = n_distinct(hgnc_symbol), .groups = "drop") %>%
    pull(n_genes)
  
  summary_list_genes <- c(
    min    = min(n_genes_vec, na.rm = TRUE),
    q1     = quantile(n_genes_vec, 0.25, na.rm = TRUE),
    median = median(n_genes_vec, na.rm = TRUE),
    mean   = mean(n_genes_vec, na.rm = TRUE),
    q3     = quantile(n_genes_vec, 0.75, na.rm = TRUE),
    max    = max(n_genes_vec, na.rm = TRUE),
    sd     = sd(n_genes_vec, na.rm = TRUE)
  )
  
  # ------------------------------
  # Final result
  # ------------------------------
  result <- list(
    dataset = dataset_name,
    summary_all  = summary_all,
    summary_up   = summary_up,
    summary_down = summary_down,
    regulation_summary = regulation_summary,
    summary_list_genes = summary_list_genes
  )
  
  class(result) <- c("gene_dataset_summary", class(result))
  message("✅ summarize_gene_dataset finished successfully for: ", dataset_name)
  return(result)
}

summarize_gene_dataset <- function(df, dataset_name = "dataset", gene_list_validation = NULL) {
  message("▶ Start summarize_gene_dataset for: ", dataset_name)
  df %>% ungroup() -> df
  stopifnot(is.data.frame(df))
  
  required_cols <- c("hgnc_symbol", "label", "source", "regulation")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0)
    stop("❌ Missing required columns: ", paste(missing_cols, collapse = ", "))
  
  # ------------------------------
  # pomocnicza: funkcja budująca tabelę walidacyjną
  # ------------------------------
  build_validation_summary <- function(data, validation_sets, mode = c("list","source")) {
    mode <- match.arg(mode)
    if (is.null(validation_sets)) return(NULL)
    
    key_cols <- if (mode == "list") c("label","hgnc_symbol") else c("source","hgnc_symbol")
    
    purrr::map_dfr(names(validation_sets), function(set_name) {
      validation_genes <- unique(validation_sets[[set_name]])
      n_total_genes <- length(validation_genes)
      
      # deduplikacja po kluczu (lista/publikacja) + gen
      genes_in_df <- data %>%
        dplyr::select(dplyr::all_of(key_cols)) %>%
        dplyr::distinct() %>%
        dplyr::filter(hgnc_symbol %in% validation_genes)
      
      if (nrow(genes_in_df) == 0) {
        return(tibble::tibble(
          validation_set = set_name,
          n_total_validation_genes = n_total_genes,
          n_unique_genes = 0,
          top_gene_1 = NA_character_, top_gene_1_count = 0,
          top_gene_2 = NA_character_, top_gene_2_count = 0,
          top_gene_3 = NA_character_, top_gene_3_count = 0,
          sum_counts = 0, mean_count = NA_real_, sd_count = NA_real_,
          fraction_detected_genes = 0
        ))
      }
      
      freq_tbl <- genes_in_df %>%
        dplyr::count(hgnc_symbol, name = "count") %>%
        dplyr::arrange(dplyr::desc(count), hgnc_symbol)
      
      top_genes  <- freq_tbl$hgnc_symbol[1:3]
      top_counts <- freq_tbl$count[1:3]
      if (length(top_genes)  < 3) top_genes  <- c(top_genes,  rep(NA, 3 - length(top_genes)))
      if (length(top_counts) < 3) top_counts <- c(top_counts, rep(0,  3 - length(top_counts)))
      
      tibble::tibble(
        validation_set = set_name,
        n_total_validation_genes = n_total_genes,
        n_unique_genes = dplyr::n_distinct(freq_tbl$hgnc_symbol),
        top_gene_1 = top_genes[1], top_gene_1_count = top_counts[1],
        top_gene_2 = top_genes[2], top_gene_2_count = top_counts[2],
        top_gene_3 = top_genes[3], top_gene_3_count = top_counts[3],
        sum_counts = sum(freq_tbl$count),
        mean_count = mean(freq_tbl$count),
        sd_count = stats::sd(freq_tbl$count),
        fraction_detected_genes = dplyr::n_distinct(freq_tbl$hgnc_symbol) / n_total_genes
      )
    }) %>%
      dplyr::arrange(dplyr::desc(top_gene_1_count))
  }
  
  # ------------------------------
  # core summarization
  # ------------------------------
  summarize_core <- function(data) {
    if (nrow(data) == 0) {
      return(list(
        n_records = 0, n_genes = 0, n_publications = 0, n_geneLists = 0,
        freq_genes_per_list  = tibble::tibble(freq = integer(), n_genes = integer(), data = list()),
        freq_genes_per_paper = tibble::tibble(freq = integer(), n_genes = integer(), data = list()),
        gene_list_validation_list   = NULL,
        gene_list_validation_source = NULL
      ))
    }
    
    result <- list(
      n_records      = nrow(data),
      n_genes        = dplyr::n_distinct(data$hgnc_symbol),
      n_publications = dplyr::n_distinct(data$source),
      n_geneLists    = dplyr::n_distinct(data$label)
    )
    
    # klasyczne podsumowania
    result$freq_genes_per_list <- data %>%
      dplyr::distinct(label, hgnc_symbol) %>%
      dplyr::count(hgnc_symbol, name = "freq") %>%
      dplyr::count(freq, name = "n_genes") %>%
      dplyr::arrange(freq)
    
    result$freq_genes_per_paper <- data %>%
      dplyr::distinct(source, hgnc_symbol) %>%
      dplyr::count(hgnc_symbol, name = "freq") %>%
      dplyr::count(freq, name = "n_genes") %>%
      dplyr::arrange(freq)
    
    # NOWE: walidacja względem zewnętrznych list genów
    result$gene_list_validation_list   <- build_validation_summary(data, gene_list_validation, mode = "list")
    result$gene_list_validation_source <- build_validation_summary(data, gene_list_validation, mode = "source")
    
    result
  }
  
  # ------------------------------
  # regulation summary
  # ------------------------------
  regulation_summary <- df %>%
    dplyr::filter(!is.na(regulation) & regulation %in% c("up", "down")) %>%
    dplyr::distinct(hgnc_symbol, regulation) %>%
    dplyr::group_by(hgnc_symbol) %>%
    dplyr::summarise(
      has_up = any(regulation == "up"),
      has_down = any(regulation == "down"),
      regulation_type = dplyr::case_when(
        has_up & has_down ~ "both",
        has_up & !has_down ~ "only_up",
        !has_up & has_down ~ "only_down",
        TRUE ~ NA_character_
      ),
      .groups = "drop"
    ) %>%
    dplyr::group_by(regulation_type) %>%
    dplyr::summarise(n_genes = dplyr::n(), genes = list(hgnc_symbol), .groups = "drop")
  
  # ------------------------------
  # global summaries
  # ------------------------------
  summary_all  <- summarize_core(df)
  summary_up   <- summarize_core(df %>% dplyr::filter(regulation == "up"))
  summary_down <- summarize_core(df %>% dplyr::filter(regulation == "down"))
  
  # ------------------------------
  # gene count distribution per list
  # ------------------------------
  n_genes_vec <- df %>%
    dplyr::group_by(label) %>%
    dplyr::summarise(n_genes = dplyr::n_distinct(hgnc_symbol), .groups = "drop") %>%
    dplyr::pull(n_genes)
  
  summary_list_genes <- c(
    min    = min(n_genes_vec, na.rm = TRUE),
    q1     = stats::quantile(n_genes_vec, 0.25, na.rm = TRUE),
    median = stats::median(n_genes_vec, na.rm = TRUE),
    mean   = base::mean(n_genes_vec, na.rm = TRUE),
    q3     = stats::quantile(n_genes_vec, 0.75, na.rm = TRUE),
    max    = max(n_genes_vec, na.rm = TRUE),
    sd     = stats::sd(n_genes_vec, na.rm = TRUE)
  )
  
  # ------------------------------
  # result
  # ------------------------------
  result <- list(
    dataset = dataset_name,
    summary_all  = summary_all,
    summary_up   = summary_up,
    summary_down = summary_down,
    regulation_summary = regulation_summary,
    summary_list_genes = summary_list_genes
  )
  
  class(result) <- c("gene_dataset_summary", class(result))
  message("✅ summarize_gene_dataset finished successfully for: ", dataset_name)
  return(result)
}

