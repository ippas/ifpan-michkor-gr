summarize_gene_dataset_and_enrichr <- function(
    df,
    dataset_name = "dataset",
    gene_list_validation = NULL,
    enrichr_run = TRUE,
    enrichr_databases = c("ChEA_2022","CellMarker_2024","LINCS_L1000_Chem_Pert_Consensus_Sigs")
) {
  message("▶ Start summarize_gene_dataset_and_enrichr for: ", dataset_name)
  df %>% dplyr::ungroup() -> df
  stopifnot(is.data.frame(df))
  
  required_cols <- c("hgnc_symbol", "label", "source", "regulation")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0)
    stop("❌ Missing required columns: ", paste(missing_cols, collapse = ", "))
  
  # -- helper: walidacja list genów
  build_validation_summary <- function(data, validation_sets, mode = c("list","source")) {
    mode <- match.arg(mode)
    if (is.null(validation_sets)) return(NULL)
    key_cols <- if (mode == "list") c("label","hgnc_symbol") else c("source","hgnc_symbol")
    
    purrr::map_dfr(names(validation_sets), function(set_name) {
      validation_genes <- unique(validation_sets[[set_name]])
      n_total_genes <- length(validation_genes)
      
      genes_in_df <- data %>%
        dplyr::select(dplyr::all_of(key_cols)) %>% dplyr::distinct() %>%
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
      
      freq_tbl <- genes_in_df %>% dplyr::count(hgnc_symbol, name = "count") %>%
        dplyr::arrange(dplyr::desc(count), hgnc_symbol)
      
      tg <- freq_tbl$hgnc_symbol[1:3]; tc <- freq_tbl$count[1:3]
      if (length(tg) < 3) tg <- c(tg, rep(NA, 3-length(tg)))
      if (length(tc) < 3) tc <- c(tc, rep(0,  3-length(tc)))
      
      tibble::tibble(
        validation_set = set_name,
        n_total_validation_genes = n_total_genes,
        n_unique_genes = dplyr::n_distinct(freq_tbl$hgnc_symbol),
        top_gene_1 = tg[1], top_gene_1_count = tc[1],
        top_gene_2 = tg[2], top_gene_2_count = tc[2],
        top_gene_3 = tg[3], top_gene_3_count = tc[3],
        sum_counts = sum(freq_tbl$count),
        mean_count = mean(freq_tbl$count),
        sd_count = stats::sd(freq_tbl$count),
        fraction_detected_genes = dplyr::n_distinct(freq_tbl$hgnc_symbol) / n_total_genes
      )
    }) %>% dplyr::arrange(dplyr::desc(top_gene_1_count))
  }
  
  # -- helper: sekcja Enrichr „w środku” summary_* (≥3/≥4/≥5 publikacji)
  run_internal_enrichr <- function(summary_obj) {
    if (!enrichr_run) return(NULL)
    fgp <- summary_obj$freq_genes_per_paper
    if (is.null(fgp) || nrow(fgp) == 0) return(NULL)
    
    genes_3pub <- fgp %>% dplyr::filter(freq >= 3) %>% dplyr::pull(data) %>% unlist() %>% unique()
    genes_4pub <- fgp %>% dplyr::filter(freq >= 4) %>% dplyr::pull(data) %>% unlist() %>% unique()
    genes_5pub <- fgp %>% dplyr::filter(freq >= 5) %>% dplyr::pull(data) %>% unlist() %>% unique()
    
    if (!requireNamespace("enrichR", quietly = TRUE)) install.packages("enrichR")
    library(enrichR)
    
    enr <- function(gset) {
      if (length(gset) == 0) return(NULL)
      enriched <- enrichR::enrichr(gset, databases = enrichr_databases)
      # filtr i top dla każdej bazy; zwracamy listę data.frameów per baza
      purrr::imap(enriched, function(tbl, db_name) {
        tbl %>% dplyr::mutate(Database = db_name) %>%
          dplyr::arrange(Adjusted.P.value) %>% 
          dplyr::mutate(
            n_genes = as.numeric(sub("/.*", "", Overlap))
          ) %>% 
          select(!c("Old.P.value", "Old.Adjusted.P.value"))
      })
    }
    
    list(
      genes_3pub = genes_3pub,
      genes_4pub = genes_4pub,
      genes_5pub = genes_5pub,
      enrichr_3pub = enr(genes_3pub),
      enrichr_4pub = enr(genes_4pub),
      enrichr_5pub = enr(genes_5pub)
    )
  }
  
  # -- rdzeń: liczenia summary_* (rozszerzone o walidację i Enrichr)
  summarize_core <- function(data) {
    if (nrow(data) == 0) {
      return(list(
        n_records = 0, n_genes = 0, n_publications = 0, n_geneLists = 0,
        freq_genes_per_list = tibble::tibble(),
        freq_genes_per_paper = tibble::tibble(),
        gene_list_validation_list = NULL,
        gene_list_validation_source = NULL,
        enrichr_summary_source = NULL
      ))
    }
    
    result <- list(
      n_records      = nrow(data),
      n_genes        = dplyr::n_distinct(data$hgnc_symbol),
      n_publications = dplyr::n_distinct(data$source),
      n_geneLists    = dplyr::n_distinct(data$label)
    )
    
    result$freq_genes_per_list <- data %>%
      dplyr::distinct(label, hgnc_symbol) %>%
      dplyr::count(hgnc_symbol, name = "freq") %>%
      dplyr::group_by(freq) %>%
      dplyr::summarise(n_genes = dplyr::n(), data = list(hgnc_symbol), .groups = "drop") %>%
      dplyr::arrange(freq)
    
    result$freq_genes_per_paper <- data %>%
      dplyr::distinct(source, hgnc_symbol) %>%
      dplyr::count(hgnc_symbol, name = "freq") %>%
      dplyr::group_by(freq) %>%
      dplyr::summarise(n_genes = dplyr::n(), data = list(hgnc_symbol), .groups = "drop") %>%
      dplyr::arrange(freq)
    
    result$gene_list_validation_list   <- build_validation_summary(data, gene_list_validation, mode = "list")
    result$gene_list_validation_source <- build_validation_summary(data, gene_list_validation, mode = "source")
    
    # << tu wstrzykujemy Enrichr DO summary_* >>
    result$enrichr_summary_source <- run_internal_enrichr(result)
    
    result
  }
  
  # -- regulation summary
  regulation_summary <- df %>%
    dplyr::filter(!is.na(regulation) & regulation %in% c("up","down")) %>%
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
  
  # -- global: all / up / down
  summary_all  <- summarize_core(df)
  summary_up   <- summarize_core(df %>% dplyr::filter(regulation == "up"))
  summary_down <- summarize_core(df %>% dplyr::filter(regulation == "down"))
  
  # -- statystyki rozmiarów list
  n_genes_vec <- df %>%
    dplyr::group_by(label) %>%
    dplyr::summarise(n_genes = dplyr::n_distinct(hgnc_symbol), .groups = "drop") %>%
    dplyr::pull(n_genes)
  
  summary_list_genes <- c(
    min = min(n_genes_vec, na.rm = TRUE),
    q1 = stats::quantile(n_genes_vec, 0.25, na.rm = TRUE),
    median = stats::median(n_genes_vec, na.rm = TRUE),
    mean = base::mean(n_genes_vec, na.rm = TRUE),
    q3 = stats::quantile(n_genes_vec, 0.75, na.rm = TRUE),
    max = max(n_genes_vec, na.rm = TRUE),
    sd = stats::sd(n_genes_vec, na.rm = TRUE)
  )
  
  result <- list(
    dataset = dataset_name,
    summary_all  = summary_all,   # ← zawiera enrichr_summary_source
    summary_up   = summary_up,    # ← zawiera enrichr_summary_source
    summary_down = summary_down,  # ← zawiera enrichr_summary_source
    regulation_summary = regulation_summary,
    summary_list_genes = summary_list_genes
  )
  
  class(result) <- c("gene_dataset_summary", class(result))
  message("✅ summarize_gene_dataset_and_enrichr finished for: ", dataset_name)
  return(result)
}
