# ##############################################################################
# ---- uses data ----
# ##############################################################################

AllBrainGeneDf
dataset_list

# ##############################################################################
dataset_list$all

dataset_list$all %>% 
  select(label, hgnc_symbol) %>% unique


read.delim(file = "results/gr-signatures/gr-signatures-multi-approach-24.04.2025.tsv") %>% 
  filter(signature_derivation == "marpiech_cluster") %>% 
  select(-signature_derivation) %>% 
  group_by(signature_name) %>% 
  summarise(genes = list(hgnc_symbol)) %>% 
  deframe() -> marpiech_clusters

map_dfr(names(marpiech_clusters), function(cluster_name) {
  cluster_genes <- unique(marpiech_clusters[[cluster_name]])     # unikalne geny w klastrze (definicja)
  n_total_cluster_genes <- length(cluster_genes)
  
  # geny z Twojego DF należące do tego klastra (do liczenia częstości)
  cluster_genes_in_df <- dataset_list$all %>%
    ungroup %>% 
    select(label, hgnc_symbol) %>% 
    unique %>% 
    filter(hgnc_symbol %in% cluster_genes)
  
  if (nrow(cluster_genes_in_df) == 0) {
    return(tibble(
      cluster = cluster_name,
      n_total_cluster_genes = n_total_cluster_genes,
      n_unique_genes = 0,
      top_gene = NA_character_,
      top_gene_count = 0
    ))
  }
  
  freq_tbl <- cluster_genes_in_df %>%
    count(hgnc_symbol, name = "count") %>%
    arrange(desc(count), hgnc_symbol)
  
  tibble(
    cluster = cluster_name,
    n_total_cluster_genes = n_total_cluster_genes,               # ⬅️ ile genów ma klaster (unikalne)
    n_unique_genes = n_distinct(cluster_genes_in_df$hgnc_symbol),# unikalne geny z klastra obecne w DF
    top_gene = freq_tbl$hgnc_symbol[1],
    top_gene_count = freq_tbl$count[1]
  )
}) %>%
  arrange(desc(top_gene_count)) %>% 
  mutate(fraction_detected_genes = n_unique_genes/n_total_cluster_genes) %>% 
  arrange(fraction_detected_genes)


map_dfr(names(marpiech_clusters), function(cluster_name) {
  cluster_genes <- unique(marpiech_clusters[[cluster_name]])
  n_total_cluster_genes <- length(cluster_genes)
  
  cluster_genes_in_df <- dataset_list$all %>%
    ungroup() %>%
    select(label, hgnc_symbol) %>%
    unique() %>%
    filter(hgnc_symbol %in% cluster_genes)
  
  if (nrow(cluster_genes_in_df) == 0) {
    return(tibble(
      cluster = cluster_name,
      n_total_cluster_genes = n_total_cluster_genes,
      n_unique_genes = 0,
      top_gene_1 = NA_character_,
      top_gene_1_count = 0,
      top_gene_2 = NA_character_,
      top_gene_2_count = 0,
      top_gene_3 = NA_character_,
      top_gene_3_count = 0
    ))
  }
  
  freq_tbl <- cluster_genes_in_df %>%
    count(hgnc_symbol, name = "count") %>%
    arrange(desc(count), hgnc_symbol)
  
  # Zabezpieczenie, jeśli klaster ma mniej niż 3 geny
  top_genes <- freq_tbl$hgnc_symbol[1:3]
  top_counts <- freq_tbl$count[1:3]
  if (length(top_genes) < 3) {
    top_genes <- c(top_genes, rep(NA, 3 - length(top_genes)))
    top_counts <- c(top_counts, rep(0, 3 - length(top_counts)))
  }
  
  tibble(
    cluster = cluster_name,
    n_total_cluster_genes = n_total_cluster_genes,
    n_unique_genes = n_distinct(cluster_genes_in_df$hgnc_symbol),
    top_gene_1 = top_genes[1],
    top_gene_1_count = top_counts[1],
    top_gene_2 = top_genes[2],
    top_gene_2_count = top_counts[2],
    top_gene_3 = top_genes[3],
    top_gene_3_count = top_counts[3]
  )
}) %>%
  mutate(fraction_detected_genes = n_unique_genes / n_total_cluster_genes) %>%
  arrange(desc(top_gene_1_count))


map_dfr(names(marpiech_clusters), function(cluster_name) {
  cluster_genes <- unique(marpiech_clusters[[cluster_name]])
  n_total_cluster_genes <- length(cluster_genes)
  
  cluster_genes_in_df <- dataset_list$all %>%
    ungroup() %>%
    select(label, hgnc_symbol) %>%
    unique() %>%
    filter(hgnc_symbol %in% cluster_genes)
  
  if (nrow(cluster_genes_in_df) == 0) {
    return(tibble(
      cluster = cluster_name,
      n_total_cluster_genes = n_total_cluster_genes,
      n_unique_genes = 0,
      top_gene_1 = NA_character_,
      top_gene_1_count = 0,
      top_gene_2 = NA_character_,
      top_gene_2_count = 0,
      top_gene_3 = NA_character_,
      top_gene_3_count = 0,
      sum_counts = 0,
      mean_count = NA_real_,
      sd_count = NA_real_
    ))
  }
  
  freq_tbl <- cluster_genes_in_df %>%
    count(hgnc_symbol, name = "count") %>%
    arrange(desc(count), hgnc_symbol)
  
  # zabezpieczenie dla klastrów z <3 genami
  top_genes <- freq_tbl$hgnc_symbol[1:3]
  top_counts <- freq_tbl$count[1:3]
  if (length(top_genes) < 3) {
    top_genes <- c(top_genes, rep(NA, 3 - length(top_genes)))
    top_counts <- c(top_counts, rep(0, 3 - length(top_counts)))
  }
  
  tibble(
    cluster = cluster_name,
    n_total_cluster_genes = n_total_cluster_genes,
    n_unique_genes = n_distinct(cluster_genes_in_df$hgnc_symbol),
    top_gene_1 = top_genes[1],
    top_gene_1_count = top_counts[1],
    top_gene_2 = top_genes[2],
    top_gene_2_count = top_counts[2],
    top_gene_3 = top_genes[3],
    top_gene_3_count = top_counts[3],
    sum_counts = sum(freq_tbl$count),
    mean_count = mean(freq_tbl$count),
    sd_count = sd(freq_tbl$count)
  )
}) %>%
  mutate(fraction_detected_genes = n_unique_genes / n_total_cluster_genes) %>%
  arrange(desc(top_gene_1_count)) %>% 
  as.data.frame() %>% 
  arrange(mean_count)
