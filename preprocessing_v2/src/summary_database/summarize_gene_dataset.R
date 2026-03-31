summarize_gene_dataset <- function(df, dataset_name = "dataset") {
  
  # ---- Core summarization function ----
  summarize_core <- function(data) {
    list(
      n_records      = nrow(data),
      n_genes        = length(unique(data$hgnc_symbol)),
      n_publications = length(unique(data$source)),
      n_geneLists    = length(unique(data$label)),
      
      freq_genes_per_list = data %>%
        select(label, hgnc_symbol) %>% 
        distinct() %>%
        pull(hgnc_symbol) %>%
        table() %>%
        as.data.frame() %>% 
        setNames(c("hgnc_symbol", "freq")) %>% 
        group_by(freq) %>% 
        nest() %>% 
        arrange(freq) %>% 
        mutate(
          n_genes = map_int(data, nrow),
          data = map(data, ~ .x %>% pull(hgnc_symbol) %>% as.character),
          fraction_genes = n_genes / sum(n_genes),
          cumulative_fraction = cumsum(fraction_genes)
        ) %>% 
        select(freq, n_genes, fraction_genes, cumulative_fraction, data),
      
      freq_genes_per_paper = data %>%
        select(source, hgnc_symbol) %>% 
        distinct() %>%
        pull(hgnc_symbol) %>%
        table() %>%
        as.data.frame() %>% 
        setNames(c("hgnc_symbol", "freq")) %>% 
        group_by(freq) %>% 
        nest() %>% 
        arrange(freq) %>% 
        mutate(
          n_genes = map_int(data, nrow),
          data = map(data, ~ .x %>% pull(hgnc_symbol) %>% as.character),
          fraction_genes = n_genes / sum(n_genes),
          cumulative_fraction = cumsum(fraction_genes)
        ) %>% 
        select(freq, n_genes, fraction_genes, cumulative_fraction, data)
    )
  }
  
  # ---- Regulation summary ----
  regulation_summary <- df %>%
    select(regulation, hgnc_symbol) %>% 
    distinct() %>% 
    group_by(hgnc_symbol) %>% 
    nest() %>% 
    mutate(n_regulation = map_int(data, nrow)) %>% 
    unnest(data) %>% 
    mutate(regulation_type = case_when(
      n_regulation == 2 ~ "both",
      n_regulation == 1 & regulation == "up" ~ "only_up",
      n_regulation == 1 & regulation == "down" ~ "only_down"
    )) %>% 
    select(-regulation, -n_regulation) %>% 
    group_by(regulation_type) %>% 
    summarise(
      n_genes = n_distinct(hgnc_symbol),
      genes = list(unique(hgnc_symbol)),
      .groups = "drop"
    )
  
  # ---- Global summaries ----
  summary_all  <- summarize_core(df)
  summary_up   <- summarize_core(df %>% filter(regulation == "up"))
  summary_down <- summarize_core(df %>% filter(regulation == "down"))
  
  # ---- Summary of gene counts per list (uses existing n_genes column + SD) ----
  n_genes_vec <- df %>%
    select(label, n_genes) %>%
    distinct() %>%
    pull(n_genes)
  
  summary_list_genes <- c(
    summary(n_genes_vec),
    "SD" = sd(n_genes_vec, na.rm = TRUE)
  )
  
  # ---- Return structured list ----
  result <- list(
    dataset = dataset_name,
    summary_all  = summary_all,
    summary_up   = summary_up,
    summary_down = summary_down,
    regulation_summary = regulation_summary,
    summary_list_genes = summary_list_genes
  )
  
  class(result) <- c("gene_dataset_summary", class(result))
  return(result)
}


summarize_gene_dataset <- function(df, dataset_name = "dataset") {
  
  # ---- Sprawdź, czy zbiór nie jest pusty ----
  if (nrow(df) == 0) {
    message("⚠️ Empty dataset: ", dataset_name)
    return(list(
      dataset = dataset_name,
      summary_all = NA,
      summary_up = NA,
      summary_down = NA,
      regulation_summary = NA,
      summary_list_genes = NA
    ))
  }
  
  # ---- Pomocnicza funkcja z zabezpieczeniem ----
  safe_freq_summary <- function(data, group_col) {
    tbl <- data %>%
      select({{ group_col }}, hgnc_symbol) %>%
      distinct() %>%
      pull(hgnc_symbol) %>%
      table() %>%
      as.data.frame()
    
    # Bezpieczna nazwa kolumn
    if (ncol(tbl) == 1) {
      colnames(tbl) <- "freq"
      tbl$hgnc_symbol <- character(0)
    } else {
      colnames(tbl) <- c("hgnc_symbol", "freq")
    }
    
    if (nrow(tbl) == 0) {
      return(tibble(
        freq = numeric(0),
        n_genes = numeric(0),
        fraction_genes = numeric(0),
        cumulative_fraction = numeric(0),
        data = list()
      ))
    }
    
    tbl %>%
      group_by(freq) %>%
      nest() %>%
      arrange(freq) %>%
      mutate(
        n_genes = map_int(data, nrow),
        data = map(data, ~ .x %>% pull(hgnc_symbol) %>% as.character),
        fraction_genes = n_genes / sum(n_genes),
        cumulative_fraction = cumsum(fraction_genes)
      ) %>%
      select(freq, n_genes, fraction_genes, cumulative_fraction, data)
  }
  
  # ---- Core summarization function ----
  summarize_core <- function(data) {
    list(
      n_records      = nrow(data),
      n_genes        = length(unique(data$hgnc_symbol)),
      n_publications = length(unique(data$source)),
      n_geneLists    = length(unique(data$label)),
      freq_genes_per_list  = safe_freq_summary(data, label),
      freq_genes_per_paper = safe_freq_summary(data, source)
    )
  }
  
  # ---- Regulation summary ----
  regulation_summary <- df %>%
    select(regulation, hgnc_symbol) %>% 
    distinct() %>% 
    group_by(hgnc_symbol) %>% 
    nest() %>% 
    mutate(n_regulation = map_int(data, nrow)) %>% 
    unnest(data) %>% 
    mutate(regulation_type = case_when(
      n_regulation == 2 ~ "both",
      n_regulation == 1 & regulation == "up" ~ "only_up",
      n_regulation == 1 & regulation == "down" ~ "only_down"
    )) %>% 
    select(-regulation, -n_regulation) %>% 
    group_by(regulation_type) %>% 
    summarise(
      n_genes = n_distinct(hgnc_symbol),
      genes = list(unique(hgnc_symbol)),
      .groups = "drop"
    )
  
  # ---- Global summaries ----
  summary_all  <- summarize_core(df)
  summary_up   <- summarize_core(df %>% filter(regulation == "up"))
  summary_down <- summarize_core(df %>% filter(regulation == "down"))
  
  # ---- Summary of gene counts per list (uses existing n_genes column + SD) ----
  n_genes_vec <- df %>%
    select(label, n_genes) %>%
    distinct() %>%
    pull(n_genes)
  
  summary_list_genes <- c(
    summary(n_genes_vec),
    "SD" = sd(n_genes_vec, na.rm = TRUE)
  )
  
  # ---- Return structured list ----
  result <- list(
    dataset = dataset_name,
    summary_all  = summary_all,
    summary_up   = summary_up,
    summary_down = summary_down,
    regulation_summary = regulation_summary,
    summary_list_genes = summary_list_genes
  )
  
  class(result) <- c("gene_dataset_summary", class(result))
  return(result)
}
