# Function to read all sheets from an Excel file into a named list
read_excel_sheets <- function(file_path) {
  # Load the readxl package
  if (!requireNamespace("readxl", quietly = TRUE)) {
    install.packages("readxl")
    library(readxl)
  }
  
  # Get the names of all sheets in the Excel file
  sheet_names <- excel_sheets(file_path)
  
  # Read each sheet into a list with names
  sheets_list <- lapply(sheet_names, function(sheet) {
    read_excel(file_path, sheet = sheet)
  })
  
  # Name the list elements with the sheet names
  names(sheets_list) <- sheet_names
  
  # Return the list of data frames
  return(sheets_list)
}


save_list_to_xlsx <- function(data_list, file_path, names_sheet = NULL) {
  # Create a new workbook
  wb <- createWorkbook()
  
  # Loop through the list and add each data frame as a new sheet
  for(i in seq_along(data_list)) {
    # Use names_sheet if provided, otherwise use names from the list or default to "Sheet_i"
    sheet_name <- if(!is.null(names_sheet)) {
      names_sheet[i]
    } else if(!is.null(names(data_list))) {
      names(data_list)[i]
    } else {
      paste0("Sheet_", i)
    }
    
    # Add a new sheet with the name
    addWorksheet(wb, sheetName = sheet_name)
    
    # Write the data frame to the sheet
    writeData(wb, sheet = sheet_name, data_list[[i]])
  }
  
  # Save the workbook to the specified file
  saveWorkbook(wb, file = file_path, overwrite = TRUE)
}

# Updated function with sorting direction
wrap_columns <- function(data, group_cols, wrap_cols, sep = "|", arrange_by = NULL, arrange_desc = FALSE) {
  if (!is.null(arrange_by)) {
    # Check if arrange_desc is a named vector matching arrange_by
    if (is.logical(arrange_desc) && length(arrange_desc) == 1) {
      if (arrange_desc) {
        data <- data %>%
          arrange(across(all_of(arrange_by), desc))
      } else {
        data <- data %>%
          arrange(across(all_of(arrange_by)))
      }
    } else if (is.logical(arrange_desc) && length(arrange_desc) == length(arrange_by)) {
      # Apply ascending or descending based on the arrange_desc vector
      sort_directions <- ifelse(arrange_desc, desc, identity)
      arrange_expressions <- map2(arrange_by, sort_directions, ~ .y(.x))
      data <- data %>%
        arrange(!!!arrange_expressions)
    }
  }
  
  data %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(n_occurrences = n(), across(all_of(wrap_cols), ~str_c(.x, collapse = sep), .names = "wrapped_{.col}")) %>%
    ungroup()
}

# Define the function with the new name
filter_and_wrap_cluster_data_by_thresholds <- function(data, p_value_threshold, fdr_threshold, top_n = NULL) {
  processed <- data %>%
    filter(p_value < p_value_threshold, fdr < fdr_threshold) %>%
    group_by(cluster, TF) %>%
    arrange(p_value)
  
  # Applying top_n if it is not NULL
  if (!is.null(top_n)) {
    processed <- processed %>% slice_min(n = top_n, order_by = p_value)
  } else {
    processed <- processed %>% slice_min(order_by = p_value, with_ties = FALSE)
  }
  
  print("good")
  
  processed %>%
    ungroup() %>%
    # select(-Term) %>%
    wrap_columns(group_cols = c("TF"), wrap_cols = c("cluster", "p_value", "fdr", "n_genes", "Genes"),
                 arrange_by = "p_value", arrange_desc = FALSE, sep = ",") %>%
    arrange(desc(n_occurences))
}


# Updated function to include explicit handling of top_n in thresholds description
filter_and_wrap_cluster_data_by_thresholds <- function(data, p_value_threshold, fdr_threshold, top_n = NULL) {
  # Create a description of the thresholds applied
  thresholds_description <- sprintf("p=%.2f;fdr=%.2f", p_value_threshold, fdr_threshold)
  if (!is.null(top_n)) {
    thresholds_description <- paste0(thresholds_description, sprintf(";top_n=%d", top_n))
  } else {
    thresholds_description <- paste0(thresholds_description, ";top_n=None")
  }
  
  processed <- data %>%
    filter(p_value < p_value_threshold, fdr < fdr_threshold) %>%
    group_by(cluster, TF) %>%
    arrange(p_value)
  
  # Applying top_n if it is not NULL
  if (!is.null(top_n)) {
    processed <- processed %>% slice_min(n = top_n, order_by = p_value)
  } else {
    processed <- processed %>% slice_min(order_by = p_value, with_ties = FALSE)
  }
  
  processed %>%
    ungroup() %>%
    wrap_columns(group_cols = c("TF"), wrap_cols = c("cluster", "p_value", "fdr", "n_genes", "Genes"),
                 arrange_by = "p_value", arrange_desc = FALSE, sep = ",") %>%
    arrange(desc(n_occurrences)) %>%
    mutate(thresholds = thresholds_description) %>%  # Add the thresholds description
    relocate(thresholds, .after = n_occurrences)  # Move the thresholds column right after n_occurrences
}

apply_thresholds_to_data <- function(data, parameters_df) {
  # Iterate over each row of the parameters dataframe
  chea_results_list <- lapply(1:nrow(parameters_df), function(i) {
    # Extract the parameters for the current iteration
    params <- parameters_df[i, ]
    
    # Group the data by 'cluster', arrange by 'p_value', and select the top n rows
    tmp_data <- data %>%
      group_by(cluster) %>%
      arrange(p_value) %>%
      slice_head(n = params[, 3]) %>%
      ungroup()
    
    print(tmp_data)
    
    # Apply the filter_and_wrap_cluster_data_by_thresholds function
    filter_and_wrap_cluster_data_by_thresholds(
      tmp_data,
      p_value_threshold = 0.05,  # Fixed threshold value
      params[, 2],               # Use parameter from the second column
      params[, 3]                # Use parameter from the third column
    )
  })
  
  # Return the results
  return(chea_results_list)
}

chea_2022_enrichr_clusters <- read_excel_sheets(file_path = "results/google-drive/enrichr/gr-dependent-transcriptional-pattern/chea_2022_filt.xlsx")

# michał
chea_2022_enrichr_clusters %>% 
  bind_rows(., .id = "cluster") %>% 
  mutate(cluster = paste0("cluster_", cluster)) %>% 
  mutate(cluster = str_replace(cluster, "cluster_UP", "UP")) %>% 
  mutate(cluster = str_replace(cluster, "cluster_DOWN", "DOWN")) %>% 
  rowwise() %>% 
  group_by(cluster, Term) %>% 
  nest() %>% 
  mutate(n_genes = map(data, ~length(convert_genes_to_vector(.x$Genes, split = ";")))) %>% 
  unnest() %>% 
  ungroup %>% 
  rename(p_value = "P-value") %>% 
  rename(fdr = "Adjusted P-value")  %>%
  select(c(cluster, Term, p_value, fdr, Genes, n_genes)) %>% 
  mutate(TF = str_split_fixed(Term, " ", 2)[, 1]) %>% 
  filter(n_genes >= 2) %>%
  filter(p_value < 0.05) %>% 
  group_by(cluster, TF) %>% 
  # filter(TF == "STAT3") %>% 
  arrange(p_value) %>% 
  slice_min(order_by = p_value, with_ties = FALSE) %>% 
  ungroup %>% 
  select(-c(Term)) -> processing_chea_TF

# mateusz
chea_df_raw %>%
  rowwise() %>% 
  group_by(cluster, Term) %>% 
  nest() %>% 
  mutate(n_genes = map(data, ~length(convert_genes_to_vector(.x$Genes, split = ";")))) %>% 
  unnest() %>% 
  ungroup %>% 
  rename(p_value = "P.value") %>% 
  rename(fdr = "Adjusted.P.value")  %>%
  select(c(cluster, Term, p_value, fdr, Genes, n_genes)) %>% 
  mutate(TF = str_split_fixed(Term, " ", 2)[, 1]) %>% 
  filter(n_genes >= 2) %>%
  filter(p_value < 0.05) %>% 
  group_by(cluster, TF) %>% 
  # filter(TF == "STAT3") %>% 
  arrange(p_value) %>% 
  slice_min(order_by = p_value, with_ties = FALSE) %>% 
  ungroup %>% 
  select(-c(Term)) -> processing_chea_TF


# processing_chea_TF %>%   
#   wrap_columns(data =., group_cols = c("TF"), wrap_cols = c("cluster", "p_value", "fdr", "n_genes", "Genes"), arrange_by = "p_value", arrange_desc = F, sep = ",") %>% 
#   arrange(desc(n_occurences))
  
filter_and_wrap_cluster_data_by_thresholds(processing_chea_TF, 0.05, 0.01)


threhold_parameters_df <- data.frame(p_value_threhold = c(0.05, 0.05, 0.05, 1, 1), fdr_threhold = c(1, 0.1, 0.01, 1, 1), top_n = c(10000, 10000, 10000, 5, 3))


# Apply the filter_and_wrap_cluster_data_by_thresholds function for each parameter set
chea_results_list <- lapply(1:nrow(threhold_parameters_df), function(i) {
  params <- threhold_parameters_df[i, ]
  
  processing_chea_TF %>%
    group_by(cluster) %>%
    arrange(p_value) %>%
    slice_head(n = params[,3]) %>%
    ungroup -> tmp_data
  
  # Run the function using the parameters from the dataframe
  filter_and_wrap_cluster_data_by_thresholds(
    # processing_chea_TF,
    tmp_data,
    p_value_threshold = params[,1],
    params[,2],
    params[,3])
})

chea_results_list %>% lapply(., dim)
  

chea_results_list %>% 
  do.call(rbind, .) %>%
  mutate(thresholds = str_replace(thresholds, "top_n=10000", "top_n=None")) %>% 
  write_tsv_xlsx(., tsv_file = "results/google-drive/enrichr/gr-dependent-transcriptional-pattern/chea2022_marpiech_clusters_summaryTF.tsv")


# chea_results_list %>%
#   lapply(., function(x) {
#     x %>% select(-thresholds)
#   }) %>%

apply_thresholds_to_data(data = processing_chea_TF, parameters_df = threhold_parameters_df ) %>% 
  lapply(., function(x) x %>% select(-thresholds)) %>% 
  save_list_to_xlsx(
    .,
    "results/google-drive/enrichr/gr-dependent-transcriptional-pattern/michkor-tables-enrichr/chea2022_marpiechcluster_summary_TF.xlsx",
    names_sheet = c("p0.05", "fdr0.1", "fdr0.01", "top5", "top3")
  )

################################################################################
# go process
go_process <- read_excel_sheets(file_path = "results/google-drive/enrichr/gr-dependent-transcriptional-pattern/michkor-tables-enrichr/GO_BP_2023_filtr.xlsx")

go_process %>% 
  discard(~ nrow(.x) == 0) %>% 
  bind_rows(., .id = "cluster") %>% 
  mutate(cluster = paste0("cluster_", cluster)) %>% 
  mutate(cluster = str_replace(cluster, "cluster_UP", "UP")) %>% 
  mutate(cluster = str_replace(cluster, "cluster_DOWN", "DOWN")) %>% 
  rowwise() %>% 
  group_by(cluster, Term) %>% 
  nest() %>% 
  mutate(n_genes = map(data, ~length(convert_genes_to_vector(.x$Genes, split = ";")))) %>% 
  unnest() %>% 
  ungroup %>% 
  rename(p_value = "P-value") %>% 
  rename(fdr = "Adjusted P-value")  %>%
  select(c(cluster, Term, p_value, fdr, Genes, n_genes)) %>% 
  # mutate(TF = str_split_fixed(Term, " ", 2)[, 1]) %>% 
  filter(n_genes >= 2) %>% 
  # filter(p_value < 0.05) %>% 
  rename(TF = "Term") -> processing_go

filter_and_wrap_cluster_data_by_thresholds(processing_go, 0.05, 0.01)

apply_thresholds_to_data(data = processing_go, parameters_df = threhold_parameters_df ) %>% 
  lapply(., function(x) {x %>% select(-thresholds) %>% rename(go_term = "TF")}) %>% 
  save_list_to_xlsx(
    .,
    "results/google-drive/enrichr/gr-dependent-transcriptional-pattern/michkor-tables-enrichr/GO_BP_2023_marpiechcluster_summary.xlsx",
    names_sheet = c("p0.05", "fdr0.1", "fdr0.01", "top5", "top3")
  )

  
################################################################################
# cellmakrer 2021

cellmarker <- read_excel_sheets(file_path = "results/google-drive/enrichr/gr-dependent-transcriptional-pattern/michkor-tables-enrichr/cellmarker_2021_filt.xlsx")

cellmarker %>% 
  discard(~ nrow(.x) == 0) %>% 
  bind_rows(., .id = "cluster") %>% 
  mutate(cluster = paste0("cluster_", cluster)) %>% 
  mutate(cluster = str_replace(cluster, "cluster_UP", "UP")) %>% 
  mutate(cluster = str_replace(cluster, "cluster_DOWN", "DOWN")) %>% 
  rowwise() %>% 
  group_by(cluster, Term) %>% 
  nest() %>% 
  mutate(n_genes = map(data, ~length(convert_genes_to_vector(.x$Genes, split = ";")))) %>% 
  unnest() %>% 
  ungroup %>% 
  rename(p_value = "P-value") %>% 
  rename(fdr = "Adjusted P-value")  %>%
  select(c(cluster, Term, p_value, fdr, Genes, n_genes)) %>% 
  # mutate(TF = str_split_fixed(Term, " ", 2)[, 1]) %>% 
  filter(n_genes >= 2) %>% 
  # filter(p_value > 0.05)
  rename(TF = "Term") -> processing_cellmarker

apply_thresholds_to_data(data = processing_cellmarker, parameters_df = threhold_parameters_df ) %>% 
  lapply(., function(x) {x %>% select(-thresholds) %>% rename(term = "TF")}) %>% 
  save_list_to_xlsx(
    .,
    "results/google-drive/enrichr/gr-dependent-transcriptional-pattern/michkor-tables-enrichr/cellmarker_marpiechcluster_summary.xlsx",
    names_sheet = c("p0.05", "fdr0.1", "fdr0.01", "top5", "top3")
  )

################################################################################
# dbsig 2021
dbsig <- read_excel_sheets(file_path = "results/google-drive/enrichr/gr-dependent-transcriptional-pattern/michkor-tables-enrichr/dbsig_filtr.xlsx")

dbsig %>% 
  discard(~ nrow(.x) == 0) %>% 
  bind_rows(., .id = "cluster") %>% 
  mutate(cluster = paste0("cluster_", cluster)) %>% 
  mutate(cluster = str_replace(cluster, "cluster_UP", "UP")) %>% 
  mutate(cluster = str_replace(cluster, "cluster_DOWN", "DOWN")) %>% 
  rowwise() %>% 
  group_by(cluster, Term) %>% 
  nest() %>% 
  mutate(n_genes = map(data, ~length(convert_genes_to_vector(.x$Genes, split = ";")))) %>% 
  unnest() %>% 
  ungroup %>% 
  rename(p_value = "P-value") %>% 
  rename(fdr = "Adjusted P-value")  %>%
  select(c(cluster, Term, p_value, fdr, Genes, n_genes)) %>% 
  # mutate(TF = str_split_fixed(Term, " ", 2)[, 1]) %>% 
  filter(n_genes >= 2) %>% 
  # filter(p_value > 0.05)
  rename(TF = "Term") -> processing_dbsig

apply_thresholds_to_data(data = processing_dbsig, parameters_df = threhold_parameters_df ) %>% 
  lapply(., function(x) {x %>% select(-thresholds) %>% rename(term = "TF")}) %>% 
  save_list_to_xlsx(
    .,
    "results/google-drive/enrichr/gr-dependent-transcriptional-pattern/michkor-tables-enrichr/dbsig_marpiechcluster_summary.xlsx",
    names_sheet = c("p0.05", "fdr0.1", "fdr0.01", "top5", "top3")
  )
