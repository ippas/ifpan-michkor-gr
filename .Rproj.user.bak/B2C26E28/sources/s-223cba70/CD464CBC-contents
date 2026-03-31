library(parallel)


enrichr_multiple_databases <- function(data_vector, databases) {
  n_databases <- length(databases)
  
  enrichr_results <- lapply(c(1:n_databases), function(i) {
    print(i)
    result <- tryCatch({
      enrichr(data_vector, databases[i])
    }, error = function(e) {
      message(paste("Error in database:", databases[i], "->", e$message))
      return(NULL)
    })
    return(result)
  })
  
  enrichr_results <- unlist(enrichr_results, recursive = FALSE)
  
  enrichr_results <- enrichr_results[sapply(enrichr_results, nrow) > 0] %>% 
    bind_rows(., .id = "database_name")
  
  return(enrichr_results)
}


enrichr_multiple_list_and_databases <- function(data_list, databases) {
  n_list <- length(data_list)
  
  # Get the names of the lists
  list_names <- names(data_list)
  
  # Apply enrichr_multiple_databases to each item in data_list
  enrichr_results <- mclapply(1:n_list, function(i) {
    print(list_names[i])  # Print the name of the current list
    enrichr_multiple_databases(data = data_list[[i]], databases = databases)
  }, mc.cores = 18)  # Parallelize across lists
  

  
  return(enrichr_results)
}

enrichr_multiple_list_and_databases <- function(data_list, databases) {
  n_list <- length(data_list)
  
  # Get the names of the lists
  list_names <- names(data_list)
  
  # Apply enrichr_multiple_databases to each item in data_list
  enrichr_results <- mclapply(1:n_list, function(i) {
    print(list_names[i])  # Print the name of the current list
    result <- enrichr_multiple_databases(data = data_list[[i]], databases = databases)
    return(list(name = list_names[i], result = result))
  }, mc.cores = 18)  # Parallelize across lists
  
  # Combine the results into a named list
  named_results <- setNames(lapply(enrichr_results, function(x) x$result), 
                            sapply(enrichr_results, function(x) x$name))
  
  return(named_results)
}

enrichr_all_database_gr_dependent_transcriptional_pattern_list <-
  enrichr_multiple_list_and_databases(data = gr_database_blocked_gene_lists$marpiech_cluster_dex,
                                      databases = dbs$libraryName)

enrichr_all_database_gr_dependent_transcriptional_pattern_list %>% 
  bind_rows(., .id = "cluster") %>% 
  # head(1000) %>%
  rename(FDR = Adjusted.P.value) %>%
  mutate(split_genes = map(Genes, ~ sort(convert_genes_to_vector(.x, split = ";")))) %>%
  mutate(Genes = map(split_genes, ~ paste(.x, collapse = ";"))) %>%
  mutate(Genes = as.character(Genes)) %>%
  select(!split_genes) %>%
  group_by(cluster, database_name, Genes) %>%
  slice_min(FDR) %>%
  slice(1) %>%
  ungroup %>%
  mutate(n_genes = as.numeric(str_split(Overlap, "/", simplify = TRUE)[, 1])) %>% 
  left_join(., n_genes_cluster, by = "cluster") -> enrichr_all_database_gr_dependent_transcriptional_pattern_df


enrichr_all_database_gr_dependent_transcriptional_pattern_df %>%
  filter(FDR < 0.05, n_genes >= 2, cluster_size <= 15) %>% as.data.frame() %>% 
  group_by(cluster, database_name) %>% 
  nest() %>% 
  mutate(mean_genes = map(data, ~mean(.x$n_genes)),
         n_term = map(data, ~nrow(.x))) %>% 
  unnest(c(mean_genes, n_term)) %>% 
  select(-data) %>% 
  as.data.frame() %>% 
  group_by(database_name) %>% 
  nest() %>% 
  mutate(mean_genes_cluster = map(data, ~mean(.x$mean_genes)),
         mean_n_term = map(data, ~mean(.x$n_term)),
         n_cluster = map(data, ~nrow(.x))) %>% 
  unnest(c(mean_genes_cluster, mean_n_term, n_cluster)) %>%
  filter(mean_n_term > 1, n_cluster > 3)



enrichr_all_database_gr_dependent_transcriptional_pattern_df %>%
  filter(FDR < 0.05, n_genes >= 2, cluster_size <= 15) %>% 
  filter(database_name == "Proteomics_Drug_Atlas_2023") %>% as.data.frame()

enrichr_all_database_gr_dependent_transcriptional_pattern_df %>%
  filter(FDR < 0.05, n_genes >= 2, cluster_size > 15) %>% 
  filter(database_name == "GWAS_Catalog_2023") %>% as.data.frame()

desired_cluster_order <- c(
  "cluster_A", "cluster_B", "cluster_C", "cluster_D", "cluster_E", 
  "cluster_F", "cluster_G", "cluster_H", "cluster_I", "cluster_J", 
  "cluster_K", "cluster_L", "cluster_M", "cluster_N", "cluster_O", 
  "cluster_P", "cluster_DOWN", "cluster_UP"
) 

enrichr_all_database_gr_dependent_transcriptional_pattern_df %>% 
  mutate(
    cluster_number = as.numeric(sub("cluster_", "", cluster)),
    new_cluster = case_when(
      cluster_number %in% 1:16 ~ paste0("cluster_", LETTERS[cluster_number]),
      cluster_number == 17 ~ "cluster_DOWN",
      cluster_number == 18 ~ "cluster_UP",
      TRUE ~ NA_character_
    ),
    cluster2 = new_cluster
  ) %>%
  mutate(cluster = factor(cluster, levels = desired_cluster_order)) -> enrichr_all_database_gr_dependent_transcriptional_pattern_df 

enrichr_all_database_gr_dependent_transcriptional_pattern_df %>% 
  select(!c(cluster_number, cluster2, new_cluster)) -> enrichr_all_database_gr_dependent_transcriptional_pattern_df

#  CellMarker_2024., ChEA_2022, GO_Biological_Process_2023
enrichr_all_database_gr_dependent_transcriptional_pattern_df %>% 
  filter(database_name %in% c("CellMarker_2024", "ChEA_2022", "GO_Biological_Process_2023")) %>% as.data.frame() %>% 
  filter(n_genes >= 2, P.value < 0.05) %>% 
  filter(database_name == "GO_Biological_Process_2023") %>%
  group_by(cluster, database_name) %>% 
  arrange(FDR, .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  as.data.frame() %>% 
  write_tsv_xlsx(., tsv_file = "results/google-drive/enrichr/gr-dependent-transcriptional-pattern/enrichr-clusters-go-bp-2023.tsv")


 

#  CellMarker_2024., ChEA_2022, GO_Biological_Process_2023
enrichr_all_database_gr_dependent_transcriptional_pattern_df %>% 
  filter(database_name %in% c("CellMarker_2024", "ChEA_2022", "GO_Biological_Process_2023")) %>% as.data.frame() %>% 
  filter(n_genes >= 2, P.value < 0.05) %>% 
  filter(database_name == "ChEA_2022") %>%
  group_by(cluster, database_name) %>% 
  arrange(FDR, .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  as.data.frame() %>% 
  write_tsv_xlsx(., tsv_file = "results/google-drive/enrichr/gr-dependent-transcriptional-pattern/enrichr-clusters-chea2022.tsv")


#  CellMarker_2024., ChEA_2022, GO_Biological_Process_2023
enrichr_all_database_gr_dependent_transcriptional_pattern_df %>% 
  filter(database_name %in% c("CellMarker_2024", "ChEA_2022", "GO_Biological_Process_2023")) %>% as.data.frame() %>% 
  filter(n_genes >= 2, P.value < 0.05) %>%
  filter(database_name == "CellMarker_2024") %>%
  group_by(cluster, database_name) %>% 
  arrange(FDR, .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  as.data.frame() %>% 
  write_tsv_xlsx(., tsv_file = "results/google-drive/enrichr/gr-dependent-transcriptional-pattern/enrichr-clusters-cellmarker2024.tsv")



enrichr_all_database_gr_dependent_transcriptional_pattern_df %>% 
  # filter(database_name %in% c("CellMarker_2024", "ChEA_2022", "GO_Biological_Process_2023")) %>% as.data.frame() %>% 
  filter(n_genes >= 2, P.value < 0.05) %>% 
  # filter(database_name == "ChEA_2022") %>%
  group_by(cluster, database_name) %>% 
  arrange(FDR, .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  as.data.frame() %>%
  write_tsv_xlsx(., tsv_file = "results/google-drive/enrichr/gr-dependent-transcriptional-pattern/enrichment-cluster-allDatabases-sortByCluster.tsv")


enrichr_all_database_gr_dependent_transcriptional_pattern_df %>% 
  # filter(database_name %in% c("CellMarker_2024", "ChEA_2022", "GO_Biological_Process_2023")) %>% as.data.frame() %>% 
  filter(n_genes >= 2, P.value < 0.05) %>% 
  # filter(database_name == "ChEA_2022") %>%
  group_by(cluster, database_name) %>% 
  arrange(FDR, .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  as.data.frame() %>% 
  arrange(database_name) %>% 
  write_tsv_xlsx(., tsv_file = "results/google-drive/enrichr/gr-dependent-transcriptional-pattern/enrichment-cluster-allDatabases-sortByDatabase.tsv")



enrichr_all_database_gr_dependent_transcriptional_pattern_df %>% 
  filter(database_name %in% c("DSigDB")) %>% as.data.frame() 
 