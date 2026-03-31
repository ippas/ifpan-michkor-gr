library(dplyr)
library(purrr)
library(tidyr)
# ---- Zdefiniuj listę datasetów do porównania ----
dataset_list <- list(
  
  # 1️⃣ Pełny zbiór
  all = AllBrainGeneDf,
  
  # 2️⃣ Listy z co najmniej 10 genami
  n10 = AllBrainGeneDf %>%
    filter(n_genes >= 10),
  
  # 3️⃣ Listy z co najmniej 50 genami
  n50 = AllBrainGeneDf %>%
    filter(n_genes >= 50),
  
  # 4️⃣ Listy z co najmniej 100 genami
  n100 = AllBrainGeneDf %>%
    filter(n_genes >= 100),
  
  # 5️⃣ Listy zawierające gen FKBP5
  FKBP5_present = AllBrainGeneDf %>%
    group_by(label) %>%
    nest() %>%
    mutate(presentFKBP5 = map_lgl(data, ~ any(.x$hgnc_symbol %in% "FKBP5"))) %>%
    filter(presentFKBP5) %>%
    unnest(data),
  
  # 6️⃣ Listy zawierające geny z cluster_P
  clusterP_present = AllBrainGeneDf %>%
    group_by(label) %>%
    nest() %>%
    mutate(presentClusterP = map_lgl(data, ~ any(.x$hgnc_symbol %in% AllBrain2BrainSignatures2GlobalMarpiechClusters$cluster_P))) %>%
    filter(presentClusterP) %>%
    unnest(data),
  
  # 7️⃣ Listy bez wartości log2ratio (brak danych ilościowych)
  missing_log2ratio = AllBrainGeneDf %>%
    filter(is.na(log2ratio)),
  
  # 8️⃣ Listy z określonym log2ratio
  with_log2ratio = AllBrainGeneDf %>%
    filter(!is.na(log2ratio)),
  
  # 9️⃣ Listy z ≥10 genami bez log2ratio
  n10_missing_log2ratio = AllBrainGeneDf %>%
    filter(n_genes >= 10) %>%
    filter(is.na(log2ratio)),
  
  # 🔟 Listy z ≥10 genami z log2ratio
  n10_with_log2ratio = AllBrainGeneDf %>%
    filter(n_genes >= 10) %>%
    filter(!is.na(log2ratio)),
  
  # 11️⃣ Listy z pomiarami po długim czasie (efekty chroniczne)
  long_time = AllBrainGeneDf %>%
    filter(time %in% c("240h", "720h", "240h_vs_720h", "672h", "3months", "3weeks", "10days")),
  
  # 12️⃣ Listy z efektami ostrymi (bez długiego czasu)
  short_time = AllBrainGeneDf %>%
    filter(!(time %in% c("240h", "720h", "240h_vs_720h", "672h", "3months", "3weeks", "10days"))),
  
  # 13️⃣ Listy z ≥10 genami o długim czasie
  n10_long_time = AllBrainGeneDf %>%
    filter(n_genes >= 10) %>%
    filter(time %in% c("240h", "720h", "240h_vs_720h", "672h", "3months", "3weeks", "10days")),
  
  # 14️⃣ Listy z ≥10 genami o krótkim czasie
  n10_short_time = AllBrainGeneDf %>%
    filter(n_genes >= 10) %>%
    filter(!(time %in% c("240h", "720h", "240h_vs_720h", "672h", "3months", "3weeks", "10days"))),
  
  # 15️⃣ Listy z FKBP5 w ostrych efektach
  FKBP5_short_time = AllBrainGeneDf %>%
    filter(!(time %in% c("240h", "720h", "240h_vs_720h", "672h", "3months", "3weeks", "10days"))) %>%
    group_by(label) %>%
    nest() %>%
    mutate(presentFKBP5 = map_lgl(data, ~ any(.x$hgnc_symbol %in% "FKBP5"))) %>%
    filter(presentFKBP5) %>%
    unnest(data),
  
  # 16️⃣ Listy z ≥10 genami i FKBP5 w ostrych efektach
  n10_FKBP5_short_time = AllBrainGeneDf %>%
    filter(n_genes >= 10) %>%
    filter(!(time %in% c("240h", "720h", "240h_vs_720h", "672h", "3months", "3weeks", "10days"))) %>%
    group_by(label) %>%
    nest() %>%
    mutate(presentFKBP5 = map_lgl(data, ~ any(.x$hgnc_symbol %in% "FKBP5"))) %>%
    filter(presentFKBP5) %>%
    unnest(data),
  
  # 17️⃣ Listy z cluster_P w ostrych efektach
  clusterP_short_time = AllBrainGeneDf %>%
    filter(!(time %in% c("240h", "720h", "240h_vs_720h", "672h", "3months", "3weeks", "10days"))) %>%
    group_by(label) %>%
    nest() %>%
    mutate(presentClusterP = map_lgl(data, ~ any(.x$hgnc_symbol %in% AllBrain2BrainSignatures2GlobalMarpiechClusters$cluster_P))) %>%
    filter(presentClusterP) %>%
    unnest(data),
  
  # 18️⃣ Listy z ≥10 genami i cluster_P w ostrych efektach
  n10_clusterP_short_time = AllBrainGeneDf %>%
    filter(n_genes >= 10) %>%
    filter(!(time %in% c("240h", "720h", "240h_vs_720h", "672h", "3months", "3weeks", "10days"))) %>%
    group_by(label) %>%
    nest() %>%
    mutate(presentClusterP = map_lgl(data, ~ any(.x$hgnc_symbol %in% AllBrain2BrainSignatures2GlobalMarpiechClusters$cluster_P))) %>%
    filter(presentClusterP) %>%
    unnest(data),
  
  # 19️⃣ Listy z cluster_P według źródła
  clusterP_by_source = AllBrainGeneDf %>%
    group_by(source) %>%
    nest() %>%
    mutate(presentClusterP = map_lgl(data, ~ any(.x$hgnc_symbol %in% AllBrain2BrainSignatures2GlobalMarpiechClusters$cluster_P))) %>%
    filter(presentClusterP) %>%
    unnest(data),
  
  # 20️⃣ Listy z cluster_P w ostrych efektach według źródła
  clusterP_short_time_by_source = AllBrainGeneDf %>%
    filter(!(time %in% c("240h", "720h", "240h_vs_720h", "672h", "3months", "3weeks", "10days"))) %>%
    group_by(source) %>%
    nest() %>%
    mutate(presentClusterP = map_lgl(data, ~ any(.x$hgnc_symbol %in% AllBrain2BrainSignatures2GlobalMarpiechClusters$cluster_P))) %>%
    filter(presentClusterP) %>%
    unnest(data),
  
  # 21️⃣ Listy z ≥10 genami, efekt ostry i cluster_P według źródła
  n10_clusterP_short_time_by_source = AllBrainGeneDf %>%
    filter(n_genes >= 10) %>%
    filter(!(time %in% c("240h", "720h", "240h_vs_720h", "672h", "3months", "3weeks", "10days"))) %>%
    group_by(source) %>%
    nest() %>%
    mutate(presentClusterP = map_lgl(data, ~ any(.x$hgnc_symbol %in% AllBrain2BrainSignatures2GlobalMarpiechClusters$cluster_P))) %>%
    filter(presentClusterP) %>%
    unnest(data),
  
  # 22️⃣ Listy z liczbą genów między 10 a 1000
  n10_1000 = AllBrainGeneDf %>%
    filter(n_genes >= 10 & n_genes <= 1000),
  
  # 23️⃣ Listy z liczbą genów między 10 a 1000 i efektami ostrymi
  n10_1000_short_time = AllBrainGeneDf %>%
    filter(n_genes >= 10 & n_genes <= 1000) %>%
    filter(!(time %in% c("240h", "720h", "240h_vs_720h", "672h", "3months", "3weeks", "10days"))),
  
  # 24️⃣ Listy z liczbą genów między 10 a 1000, efektami ostrymi i zawierające FKBP5 (grupowane po publikacji)
  n10_1000_short_time_FKBP5_by_source = AllBrainGeneDf %>%
    filter(n_genes >= 10 & n_genes <= 1000) %>%
    filter(!(time %in% c("240h", "720h", "240h_vs_720h", "672h", "3months", "3weeks", "10days"))) %>%
    group_by(source) %>%
    nest() %>%
    mutate(has_FKBP5 = map_lgl(data, ~ any(.x$hgnc_symbol %in% "FKBP5"))) %>%
    filter(has_FKBP5) %>%
    unnest(data),
  
  # 25️⃣ Listy z liczbą genów między 10 a 1000, efektami ostrymi i zawierające ≥1 gen z cluster_P (grupowane po publikacji)
  n10_1000_short_time_clusterP_by_source = AllBrainGeneDf %>%
    filter(n_genes >= 10 & n_genes <= 1000) %>%
    filter(!(time %in% c("240h", "720h", "240h_vs_720h", "672h", "3months", "3weeks", "10days"))) %>%
    group_by(source) %>%
    nest() %>%
    mutate(has_clusterP = map_lgl(data, ~ any(.x$hgnc_symbol %in% AllBrain2BrainSignatures2GlobalMarpiechClusters$cluster_P))) %>%
    filter(has_clusterP) %>%
    unnest(data),
  
  # 26️⃣ Listy z liczbą genów między 10 a 1000, efektami ostrymi i zawierające FKBP5 oraz ≥1 gen z cluster_P (grupowane po publikacji)
  n10_1000_short_time_clusterP_by_source = AllBrainGeneDf %>%
    filter(n_genes >= 10 & n_genes <= 1000) %>%
    filter(!(time %in% c("240h", "720h", "240h_vs_720h", "672h", "3months", "3weeks", "10days"))) %>%
    group_by(source) %>%
    nest() %>%
    mutate(has_clusterP = map_lgl(data, ~ any(.x$hgnc_symbol %in% AllBrain2BrainSignatures2GlobalMarpiechClusters$cluster_P))) %>%
    filter(has_clusterP) %>%
    unnest(data)
)


# ---- Uruchomienie funkcji dla każdego datasetu ----
gene_summary_list <- imap(dataset_list, ~ {
  message("🔍 Processing dataset: ", .y)
  summarize_gene_dataset(df = .x, dataset_name = .y)
})

gene_summary_list$all$summary_list_genes

rm(dataset_list)
summary_table <- imap_dfr(gene_summary_list, function(x, name) {
  
  # jeśli dataset jest pusty, pomiń
  if (is.null(x) || is.na(x) || all(is.na(x))) {
    message("⚠️ Skipping empty dataset: ", name)
    return(NULL)
  }
  
  # ---- summary_all ----
  n_records      <- x$summary_all$n_records
  n_genes        <- x$summary_all$n_genes
  n_publications <- x$summary_all$n_publications
  n_geneLists    <- x$summary_all$n_geneLists
  
  freq_tbl <- x$summary_all$freq_genes_per_list
  n_genes_freq1 <- if (!is.null(freq_tbl) && nrow(freq_tbl) > 0) freq_tbl %>% filter(freq == 1) %>% pull(n_genes) %>% first() else NA
  last_row <- if (!is.null(freq_tbl) && nrow(freq_tbl) > 0) freq_tbl %>% arrange(freq) %>% tail(1) else tibble(freq = NA, n_genes = NA, data = list(NA))
  max_freq <- if (!is.null(last_row$freq)) first(last_row$freq) else NA
  n_genes_at_max_freq <- if (!is.null(last_row$n_genes)) first(last_row$n_genes) else NA
  genes_at_max_freq <- if (!is.null(last_row$data[[1]])) last_row$data[[1]] %>% head(3) %>% paste(collapse = ", ") else NA
  
  # ---- regulation_summary ----
  reg_tbl <- x$regulation_summary
  both_n      <- if (!is.null(reg_tbl)) reg_tbl %>% filter(regulation_type == "both") %>% pull(n_genes) %>% first() else NA
  only_up_n   <- if (!is.null(reg_tbl)) reg_tbl %>% filter(regulation_type == "only_up") %>% pull(n_genes) %>% first() else NA
  only_down_n <- if (!is.null(reg_tbl)) reg_tbl %>% filter(regulation_type == "only_down") %>% pull(n_genes) %>% first() else NA
  
  # ---- summary_up ----
  n_records_up      <- x$summary_up$n_records
  n_genes_up        <- x$summary_up$n_genes
  n_publications_up <- x$summary_up$n_publications
  n_geneLists_up    <- x$summary_up$n_geneLists
  
  freq_tbl_up <- x$summary_up$freq_genes_per_list
  n_genes_freq1_up <- if (!is.null(freq_tbl_up) && nrow(freq_tbl_up) > 0) freq_tbl_up %>% filter(freq == 1) %>% pull(n_genes) %>% first() else NA
  last_row_up <- if (!is.null(freq_tbl_up) && nrow(freq_tbl_up) > 0) freq_tbl_up %>% arrange(freq) %>% tail(1) else tibble(freq = NA, n_genes = NA, data = list(NA))
  max_freq_up <- if (!is.null(last_row_up$freq)) first(last_row_up$freq) else NA
  n_genes_at_max_freq_up <- if (!is.null(last_row_up$n_genes)) first(last_row_up$n_genes) else NA
  genes_at_max_freq_up <- if (!is.null(last_row_up$data[[1]])) last_row_up$data[[1]] %>% head(3) %>% paste(collapse = ", ") else NA
  
  # ---- summary_down ----
  n_records_down      <- x$summary_down$n_records
  n_genes_down        <- x$summary_down$n_genes
  n_publications_down <- x$summary_down$n_publications
  n_geneLists_down    <- x$summary_down$n_geneLists
  
  freq_tbl_down <- x$summary_down$freq_genes_per_list
  n_genes_freq1_down <- if (!is.null(freq_tbl_down) && nrow(freq_tbl_down) > 0) freq_tbl_down %>% filter(freq == 1) %>% pull(n_genes) %>% first() else NA
  last_row_down <- if (!is.null(freq_tbl_down) && nrow(freq_tbl_down) > 0) freq_tbl_down %>% arrange(freq) %>% tail(1) else tibble(freq = NA, n_genes = NA, data = list(NA))
  max_freq_down <- if (!is.null(last_row_down$freq)) first(last_row_down$freq) else NA
  n_genes_at_max_freq_down <- if (!is.null(last_row_down$n_genes)) first(last_row_down$n_genes) else NA
  genes_at_max_freq_down <- if (!is.null(last_row_down$data[[1]])) last_row_down$data[[1]] %>% head(3) %>% paste(collapse = ", ") else NA
  
  # ---- summary_list_genes ----
  slg <- x$summary_list_genes
  min_n_genes_list    <- if (!is.null(slg)) slg["Min."] else NA
  q1_n_genes_list     <- if (!is.null(slg)) slg["1st Qu."] else NA
  median_n_genes_list <- if (!is.null(slg)) slg["Median"] else NA
  mean_n_genes_list   <- if (!is.null(slg)) slg["Mean"] else NA
  q3_n_genes_list     <- if (!is.null(slg)) slg["3rd Qu."] else NA
  max_n_genes_list    <- if (!is.null(slg)) slg["Max."] else NA
  sd_n_genes_list     <- if (!is.null(slg)) slg["SD"] else NA
  
  # ---- Final tibble ----
  tibble(
    dataset = name,
    n_records, n_genes, n_publications, n_geneLists,
    n_genes_freq1, max_freq, n_genes_at_max_freq, genes_at_max_freq,
    both_n, only_up_n, only_down_n,
    n_records_up, n_genes_up, n_publications_up, n_geneLists_up,
    n_genes_freq1_up, max_freq_up, n_genes_at_max_freq_up, genes_at_max_freq_up,
    n_records_down, n_genes_down, n_publications_down, n_geneLists_down,
    n_genes_freq1_down, max_freq_down, n_genes_at_max_freq_down, genes_at_max_freq_down,
    min_n_genes_list, q1_n_genes_list, median_n_genes_list, mean_n_genes_list,
    q3_n_genes_list, max_n_genes_list, sd_n_genes_list
  )
})


summary_table %>% 
  arrange(desc(n_records)) %>% 
  as.data.frame() %>% 
  mutate(percent_both = both_n/n_genes *100)
