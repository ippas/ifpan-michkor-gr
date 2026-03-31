AllBrainGeneDf %>% 
  filter(hgnc_symbol %in% AllBrain2BrainSignatures2GlobalMarpiechClusters$cluster_P) %>% 
  select(label, hgnc_symbol) %>% unique %>% 
  .$label %>% table() %>% as.data.frame %>% 
  filter(Freq > 5)

AllBrain2BrainSignatures2GlobalMarpiechClusters$cluster_M %>% length()
AllBrain2BrainSignatures2GlobalMarpiechClusters$cluster_N %>% length()
AllBrain2BrainSignatures2GlobalMarpiechClusters$cluster_P %>% length()

cluster_names <- c(
  paste0("cluster_", LETTERS[1:16]),
  paste0("cluster_", LETTERS[1:16], "_UP"),
  paste0("cluster_", LETTERS[1:16], "_DOWN")
)


map_dfr(cluster_names, function(cluster_name) {
  if (!cluster_name %in% names(AllBrain2BrainSignatures2GlobalMarpiechClusters)) return(NULL)
  
  cluster_genes <- AllBrain2BrainSignatures2GlobalMarpiechClusters[[cluster_name]]
  if (length(cluster_genes) == 0) return(NULL)
  
  df_filtered <- dataset_list$all %>%
    filter(n_genes >= 10) %>%
    filter(!(time %in% c("240h","720h","240h_vs_720h","240_vs_720h",
                         "672h","3months","3weeks","10days"))) %>%
    select(label, hgnc_symbol) %>%
    filter(hgnc_symbol %in% cluster_genes) %>%
    distinct()
  
  if (nrow(df_filtered) == 0) return(NULL)
  
  freq_table <- sort(table(df_filtered$hgnc_symbol), decreasing = TRUE)
  freq_values <- as.numeric(freq_table)
  top_genes <- names(freq_table)[1:3] %>% replace_na("")  # top 3
  
  tibble(
    cluster = cluster_name,
    direction = case_when(
      grepl("_UP$", cluster_name) ~ "UP",
      grepl("_DOWN$", cluster_name) ~ "DOWN",
      TRUE ~ "base"
    ),
    n_genes_cluster  = length(cluster_genes),
    n_detected_genes = length(freq_values),
    median_raw = median(freq_values),
    mean_raw   = mean(freq_values),
    max_raw    = max(freq_values),
    top_gene_1 = top_genes[1],
    top_gene_2 = top_genes[2],
    top_gene_3 = top_genes[3]
  )
}) %>%
  arrange(cluster, desc(median_raw))


map_dfr(cluster_names, function(cluster_name) {
  if (!cluster_name %in% names(AllBrain2BrainSignatures2GlobalMarpiechClusters)) return(NULL)
  
  cluster_genes <- AllBrain2BrainSignatures2GlobalMarpiechClusters[[cluster_name]]
  if (length(cluster_genes) == 0) return(NULL)
  
  df_filtered <- dataset_list$n10_short_time %>%
    filter(n_genes >= 10) %>%
    filter(!(time %in% c("240h","720h","240h_vs_720h","240_vs_720h",
                         "672h","3months","3weeks","10days"))) %>%
    select(label, hgnc_symbol) %>%
    filter(hgnc_symbol %in% cluster_genes) %>%
    distinct()
  
  if (nrow(df_filtered) == 0) return(NULL)
  
  freq_table <- sort(table(df_filtered$hgnc_symbol), decreasing = TRUE)
  freq_values <- as.numeric(freq_table)
  top_genes <- names(freq_table)[1:3] %>% replace_na("")  # top 3
  
  tibble(
    cluster = cluster_name,
    direction = case_when(
      grepl("_UP$", cluster_name) ~ "UP",
      grepl("_DOWN$", cluster_name) ~ "DOWN",
      TRUE ~ "base"
    ),
    n_genes_cluster  = length(cluster_genes),
    n_detected_genes = length(freq_values),
    median_raw = median(freq_values),
    mean_raw   = mean(freq_values),
    max_raw    = max(freq_values),
    top_gene_1 = top_genes[1],
    top_gene_2 = top_genes[2],
    top_gene_3 = top_genes[3]
  )
}) %>%
  arrange(cluster, desc(median_raw)) %>% 
  mutate(prop_detected_genes = n_detected_genes/n_genes_cluster) 


dataset_list$all %>%
  filter(n_genes >= 10) %>%
  filter(!(time %in% c("240h","720h","240h_vs_720h","240_vs_720h",
                       "672h","3months","3weeks","10days"))) %>% filter(hgnc_symbol == "PEX11A")




AllBrainGeneDf %>% filter(hgnc_symbol == "NTS")
  
dataset_list$n10_short_time

gene_summary_list$n10_short_time$summary_list_genes

gene_summary_list$n10_short_time$summary_up$freq_genes_per_paper  %>% filter(freq >= 4) %>% .$data %>% unlist %>% rev %>% cat(sep = "\n")

gene_summary_list$n10_short_time$summary_down$freq_genes_per_paper  %>% filter(freq >= 4) %>% .$data %>% unlist %>% rev %>% cat(sep = "\n")


intersect(
  gene_summary_list$n10_short_time$summary_up$freq_genes_per_paper %>%
    filter(freq >= 3) %>%
    .$data %>%
    unlist(),
  gene_summary_list$n10_short_time$summary_down$freq_genes_per_paper %>%
    filter(freq >= 3) %>%
    .$data %>%
    unlist()
) %>%
  cat(sep = "\n")
# super, brak wspólnych genów, super filt na przynajmniej w trzech publikacjach 