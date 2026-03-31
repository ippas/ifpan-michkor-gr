cluster_mapper <- c(
  "cluster_1" = "cluster_A",
  "cluster_2" = "cluster_B",
  "cluster_3" = "cluster_C",
  "cluster_4" = "cluster_D",
  "cluster_5" = "cluster_E",
  "cluster_6" = "cluster_F",
  "cluster_7" = "cluster_G",
  "cluster_8" = "cluster_H",
  "cluster_9" = "cluster_I",
  "cluster_10" = "cluster_J",
  "cluster_11" = "cluster_K",
  "cluster_12" = "cluster_L",
  "cluster_13" = "cluster_M",
  "cluster_14" = "cluster_N",
  "cluster_15" = "cluster_O",
  "cluster_16" = "cluster_P",
  "cluster_17" = "cluster_DOWN",
  "cluster_18" = "cluster_UP"
)

rename_list_elements <- function(list, mapper) {
  names(list) <- mapper[names(list)]
  return(list)
}



gr_database_blocked_gene_lists$marpiech_cluster_dex_letters <- rename_list_elements(gr_database_blocked_gene_lists$marpiech_cluster_dex, cluster_mapper)

gr_database_blocked_gene_lists$marpiech_cluster_dex_letters$cluster_P
gr_database_blocked_gene_lists$marpiech_cluster_dex$cluster_16

save_list_to_excel <- function(data_list, file_name) {
  # Load the necessary package
  if (!require(openxlsx)) {
    install.packages("openxlsx")
    library(openxlsx)
  }
  
  # Create a new workbook
  wb <- createWorkbook()
  
  # Loop through each element in the list and add it to the workbook
  for (i in seq_along(data_list)) {
    sheet_name <- names(data_list)[i]
    if (is.null(sheet_name) || sheet_name == "") {
      sheet_name <- paste("Sheet", i)
    }
    addWorksheet(wb, sheetName = sheet_name)
    writeData(wb, sheet = sheet_name, data_list[[i]])
  }
  
  # Save the workbook to a file
  saveWorkbook(wb, file_name, overwrite = TRUE)
  
  # Return a message confirming the file has been saved
  return(paste("Workbook saved as", file_name))
}

save_list_to_excel <- function(data_list, file_name) {
  # Load the necessary package
  if (!require(openxlsx)) {
    install.packages("openxlsx")
    library(openxlsx)
  }
  
  # Define the order of the cluster sheets
  cluster_order <- c("cluster_A", "cluster_B", "cluster_C", "cluster_D", "cluster_E", 
                     "cluster_F", "cluster_G", "cluster_H", "cluster_I", "cluster_J", 
                     "cluster_K", "cluster_L", "cluster_M", "cluster_N", "cluster_O", 
                     "cluster_P", "cluster_DOWN", "cluster_UP")
  
  # Check for and remove any missing entries from the order to avoid errors
  valid_clusters <- cluster_order %in% names(data_list)
  if (any(!valid_clusters)) {
    warning("Some specified clusters do not exist in the data list and will be skipped: ", 
            paste(cluster_order[!valid_clusters], collapse = ", "))
    cluster_order <- cluster_order[valid_clusters]
  }
  
  # Sort the list according to the predefined order
  data_list <- data_list[cluster_order]
  
  # Create a new workbook
  wb <- createWorkbook()
  
  # Loop through each element in the sorted list and add it to the workbook
  for (i in seq_along(data_list)) {
    sheet_name <- names(data_list)[i]
    # Ensure that a valid sheet name is used
    if (is.null(sheet_name) || sheet_name == "") {
      sheet_name <- paste("Sheet", i)
    }
    addWorksheet(wb, sheetName = sheet_name)
    writeData(wb, sheet = sheet_name, data_list[[i]])
  }
  
  # Save the workbook to a file
  saveWorkbook(wb, file_name, overwrite = TRUE)
  
  # Return a message confirming the file has been saved
  return(paste("Workbook saved as", file_name))
}


gr_database_blocked_gene_lists$marpiech_cluster_dex

gr_database_blocked_gene_lists$marpiech_cluster_dex_letters <- gr_database_blocked_gene_lists$marpiech_cluster_dex


  
enrichr_databases <-
  c(
    "BioPlanet_2019",
    "CellMarker_2024",
    "GO_Biological_Process_2023"
    # "KEGG_2021_Human",
    # "WikiPathway_2023_Human",
    # "Elsevier_Pathway_Collection"
  )


dbs %>% 
  filter(libraryName %in% c(enrichr_databases)) -> dbs_filtered

dbs_filtered %>% dim


enrichr_multiple_databases <- function(data_vector, databases){
  n_databases <- length(databases)
  
  lapply(c(1:n_databases), function(i){
    print(i)
    enrichr(data_vector, databases[i])
  }) -> enrichr_results
  
  
  # enrichr_results <- unlist(enrichr_results, recursive = FALSE)
  
  enrichr_results <- unlist(enrichr_results, recursive = FALSE) %>%
    bind_rows(., .id = "enrichr_database")
  
  return(enrichr_results)
}

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
  
  # Uncomment this line if you want to bind rows of the results into a single data frame
  # enrichr_results <- unlist(enrichr_results, recursive = FALSE) %>%
  #   bind_rows(., .id = "enrichr_database")
  
  return(enrichr_results)
}



enrichr_multiple_list_and_databases <- function(data_list, databases) {
  n_list <- length(data_list)
  
  # Get the names of the lists
  list_names <- names(data_list)
  
  # Apply enrichr_multiple_databases to each item in data_list
  enrichr_results <- lapply(1:n_list, function(i) {
    print(list_names[i])  # Print the name of the current list
    enrichr_multiple_databases(data = data_list[[i]],  # Extract the data, not a sublist
                               databases = databases)
  })
  
  # Assign names to the results
  names(enrichr_results) <- list_names
  
  
  return(enrichr_results)
}

all_ruslt_enrichr_list2 = enrichr_multiple_list_and_databases(data_list = gr_database_blocked_gene_lists$marpiech_cluster_dex_letters, databases = dbs$libraryName)

enrichr_all_database_gr_dependent_transcriptional_pattern_list <-
  enrichr_multiple_list_and_databases(data = gr_database_blocked_gene_lists$marpiech_cluster_dex,
                                      databases = dbs$libraryName)


enrichr_all_database_gr_dependent_transcriptional_pattern_list %>% 
  bind_rows(., .id = "cluster") %>% 
  mutate(cluster = case_when(
    cluster == "cluster_1" ~ "cluster_A",
    cluster == "cluster_2" ~ "cluster_B",
    cluster == "cluster_3" ~ "cluster_C",
    cluster == "cluster_4" ~ "cluster_D",
    cluster == "cluster_5" ~ "cluster_E",
    cluster == "cluster_6" ~ "cluster_F",
    cluster == "cluster_7" ~ "cluster_G",
    cluster == "cluster_8" ~ "cluster_H",
    cluster == "cluster_9" ~ "cluster_I",
    cluster == "cluster_10" ~ "cluster_J",
    cluster == "cluster_11" ~ "cluster_K",
    cluster == "cluster_12" ~ "cluster_L",
    cluster == "cluster_13" ~ "cluster_M",
    cluster == "cluster_14" ~ "cluster_N",
    cluster == "cluster_15" ~ "cluster_O",
    cluster == "cluster_16" ~ "cluster_P",
    cluster == "cluster_17" ~ "cluster_DOWN",
    cluster == "cluster_18" ~ "cluster_UP",
  )) ->  enrichr_all_database_gr_dependent_transcriptional_pattern_df


enrichr_all_database_gr_dependent_transcriptional_pattern_list$cluster_2 %>% filter(database_name == "DSigDB")

write_tsv_xlsx(enrichr_all_database_gr_dependent_transcriptional_pattern_df,
               tsv_file = "results/google-drive/enrichr/gr-dependent-transcriptional-pattern/enrichr-all-database-clusters-no-filters.tsv")
  
enrichr_all_database_gr_dependent_transcriptional_pattern_df %>% 
  filter(database_name == "DSigDB") %>% 
  .$cluster %>% unique()

enrichr_all_database_gr_dependent_transcriptional_pattern_df %>% .$cluster %>% unique

tmp <- Filter(function(x) nrow(x) == 0, enrichr_all_database_gr_dependent_transcriptional_pattern_list)

dsigndb_list <-
  enrichr_multiple_list_and_databases(data = gr_database_blocked_gene_lists$marpiech_cluster_dex,
                                      databases = c("DSigDB"))


dsigndb_list %>% 
  bind_rows(., .id = "cluster") %>% 
  mutate(cluster = case_when(
    cluster == "cluster_1" ~ "cluster_A",
    cluster == "cluster_2" ~ "cluster_B",
    cluster == "cluster_3" ~ "cluster_C",
    cluster == "cluster_4" ~ "cluster_D",
    cluster == "cluster_5" ~ "cluster_E",
    cluster == "cluster_6" ~ "cluster_F",
    cluster == "cluster_7" ~ "cluster_G",
    cluster == "cluster_8" ~ "cluster_H",
    cluster == "cluster_9" ~ "cluster_I",
    cluster == "cluster_10" ~ "cluster_J",
    cluster == "cluster_11" ~ "cluster_K",
    cluster == "cluster_12" ~ "cluster_L",
    cluster == "cluster_13" ~ "cluster_M",
    cluster == "cluster_14" ~ "cluster_N",
    cluster == "cluster_15" ~ "cluster_O",
    cluster == "cluster_16" ~ "cluster_P",
    cluster == "cluster_17" ~ "cluster_DOWN",
    cluster == "cluster_18" ~ "cluster_UP",
  )) -> dsigndb_df
  
new_dsigndb_list <- split(dsigndb_df, dsigndb_df$cluster) %>% 
  lapply(., function(x){x %>% select(-cluster)})

save_list_to_excel(data_list = new_dsigndb_list, file_name = "results/google-drive/enrichr/gr-dependent-transcriptional-pattern/enrichr-dsigdb-cluster-no-filters.xlsx")

data.frame(do.call(cbind, gr_database_blocked_gene_lists$marpiech_cluster_dex_letters)) %>% 
  select(unname(cluster_mapper)) 
group_by(label) %>%
  mutate(row_id = row_number()) %>% 
  ungroup() 

gr_database_blocked_gene_lists$marpiech_cluster_dex_letters %>%
  lapply(., function(x){x %>% as.data.frame() %>% set_colnames("hgnc_symbol")}) %>% 
  bind_rows(., .id="cluster") %>% 
  group_by(cluster) %>% 
  mutate(row_id = row_number()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = cluster, values_from = hgnc_symbol)  %>% 
  select(-row_id) %>% 
  select(unname(cluster_mapper)) %>% 
  write_tsv_xlsx(tsv_file = "results/google-drive/gene-list-datasets/specific-profiles-gr.tsv")


