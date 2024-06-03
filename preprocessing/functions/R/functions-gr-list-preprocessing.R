# Define the function with an additional `tissue_column` argument
categorize_tissue <- function(df, tissue_column) {
  # Check if the specified tissue_column exists in the dataframe
  if(!tissue_column %in% names(df)) {
    stop(paste("The dataframe does not contain a", tissue_column, "column."), call. = FALSE)
  }
  
  # Perform the mutations using the specified column for tissue information
  df <- df %>%
    # mutate(system = case_when(
    #   str_detect(!!sym(tissue_column), "brain|cortex(?<!adrenal-cortex)|hippocampal|neuronal|astrocyte|striatum|hypothalamus|hippocampus|ARC|pituitary|astrocyte") ~ "brain_and_nervous_system",
    #   str_detect(!!sym(tissue_column), "blood|placenta") ~ "blood_and_immune_system",
    #   str_detect(!!sym(tissue_column), "lung") ~ "respiratory_system",
    #   str_detect(!!sym(tissue_column), "liver|kidney|small-intestine|spleen") ~ "digestive_system",
    #   str_detect(!!sym(tissue_column), "muscle|adipose|adipocyte|bone|cartilage|anterior") ~ "musculoskeletal_system",
    #   str_detect(!!sym(tissue_column), "embryos") ~ "embryos",
    #   str_detect(!!sym(tissue_column), "adrenal-gland|adrenal-cortex") ~ "endocrine_system",
    #   TRUE ~ "other_tissues"
    # )) %>%
    mutate(simple_tissue = case_when(
      str_detect(!!sym(tissue_column), "brain|cortex(?<!adrenal-cortex)|hippocampal|neuronal|astrocyte|striatum|hypothalamus|hippocampus|ARC|pituitary|astrocyte") ~ "brain",
      grepl("embryos", !!sym(tissue_column), ignore.case = TRUE) ~ "embryos",
      grepl("placenta", !!sym(tissue_column), ignore.case = TRUE) ~ "placenta",
      grepl("bone", !!sym(tissue_column), ignore.case = TRUE) ~ "bone",
      grepl("cartilage", !!sym(tissue_column), ignore.case = TRUE) ~ "cartilage",
      grepl("liver", !!sym(tissue_column), ignore.case = TRUE) ~ "liver",
      grepl("kidney", !!sym(tissue_column), ignore.case = TRUE) ~ "kidney",
      grepl("spleen", !!sym(tissue_column), ignore.case = TRUE) ~ "spleen",
      grepl("adipose|adipocyte", !!sym(tissue_column), ignore.case = TRUE) ~ "adipose",
      grepl("blood", !!sym(tissue_column), ignore.case = TRUE) ~ "blood",
      grepl("lung", !!sym(tissue_column), ignore.case = TRUE) ~ "lung",
      grepl("adrenal-gland|adrenal-cortex", !!sym(tissue_column), ignore.case = TRUE) ~ "adrenal-gland",
      grepl("small-intestine", !!sym(tissue_column), ignore.case = TRUE) ~ "small-intestine",
      grepl("muscle|anterior", !!sym(tissue_column), ignore.case = TRUE) ~ "muscle",
      TRUE ~ "Other"  # This will categorize any tissue not matched above as "Other"
    )) %>% 
    mutate(detailed_tissue = case_when(
      simple_tissue == "lung" & (cell == "NA" | is.na(cell)) ~ "lung",
      simple_tissue == "lung" & grepl("ASM", cell, ignore.case = TRUE) ~ "lung-AirwaySmoothMuscle",
      simple_tissue == "lung" & grepl("A549|BEAS-2B|pHBECs|H1944|H1975|H2122|H460", cell, ignore.case = TRUE) ~ "lung-epithelial",
      simple_tissue == "blood" & grepl("THP-1|macrophages|mBMDM|hMDM", cell, ignore.case = TRUE) ~ "blood-macrophages",
      simple_tissue == "blood" & grepl("Tcell|NKcell|Monocyte|Bcell|REH-overexpression-GCR|NALM6", cell, ignore.case = TRUE) ~ "blood-noMacrophages",
      simple_tissue == "blood" & (cell == "NA" | is.na(cell)) ~ "blood",
      tissue == "prefrontal-cortex" ~ "brain-cortex",
      tissue == "ARC" ~ "brain-hypothalamus",
      tissue == "hippocampus" ~ "brain-hippocampus",
      tissue == "hippocampla-slices" ~ "brain-hippocampus",
      tissue == "hippocampal-progenitor-cell-line" ~ "brain-hippocampus",
      tissue == "hippocampal-slices" ~ "brain-hippocampus",
      tissue == "cells-derived-from-human-fetal-brain-tissue" ~ "brain-hippocampus",
      tissue == "cortex" ~ "brain-cortex",
      tissue == "neuronal-pc12-cells" ~ "brain-neurons",
      tissue == "hypothalamus" ~ "brain-hypothalamus",
      grepl("neurones|neurons|neuron", cell, ignore.case = TRUE) ~ "brain-neurons", 
      grepl("mglia|astrocytes|astrocyte|oligodendrocytes|OPC", cell, ignore.case = TRUE) ~ "brain-glia",
      TRUE ~ simple_tissue
    )) 
  
  return(df)
}


################################################################################
process_tissue_data <- function(data, columns, tissue_column = "tissue", additional_columns = NULL) {
  library(dplyr)
  
  data_processed <- data %>%
    mutate(!!tissue_column := ifelse(is.na(!!sym(tissue_column)), cell, !!sym(tissue_column))) %>%
    categorize_tissue(tissue_column = tissue_column) %>%
    select(c(columns, additional_columns)) %>%
    unique()
  
  return(data_processed)
}

barplot_tissue_data <- function(data, count_column = "simple_tissue", threshold = 3) {
  library(ggplot2)
  
  data %>%
    .[[count_column]] %>%
    table() %>%
    as.data.frame() %>% 
    set_colnames(c("tissue", "freq")) %>% 
    arrange(desc(freq)) %>%
    ggplot(aes(x = reorder(tissue, -freq), y = freq)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = threshold, color = "red", linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

################################################################################
# filter_gene_list <- function(list, ..., log2ratio = FALSE){
#   list$metadata %>% 
#     filter(...) -> filter_metadata
#   
#   filter_metadata$label -> label_vector
#   
#   list$gene_list[label_vector] -> gene_lists
#   
#   if(log2ratio == TRUE){
#     list$log2ratio[label_vector] -> log2ratio_lists
#     
#     
#     results <- list(metadata = filter_metadata,
#                     gene_lists = gene_lists,
#                     log2ratio_lists = log2ratio_lists)
#     
#   } else {
#     results <- list(metadata = filter_metadata,
#                     gene_lists = gene_lists)
#   }
#   
#   return(results)
# }

filter_gene_list <- function(list, ...) {
  # Filter metadata and create a label vector
  filter_metadata <- dplyr::filter(list$metadata, ...)
  label_vector <- filter_metadata$label
  
  # Initialize results with filtered metadata
  results <- list(metadata = filter_metadata)
  
  # Loop over all elements in the list except 'metadata'
  for (name in setdiff(names(list), "metadata")) {
    # Apply the filter based on the label_vector if the element is a vector or a matrix
    if (is.vector(list[[name]]) || is.matrix(list[[name]])) {
      results[[name]] <- list[[name]][label_vector]
    }
  }
  
  return(results)
}

################################################################################
convert_genes_to_vector <- function(data, split = "\\|"){
  lapply(data, strsplit, split = split) %>% unlist
}

################################################################################
# Function to process gene lists and return results in a table
calculate_gene_lists_diversity_metrics <- function(gene_lists) {
  # Calculating Gini Coefficient
  gini_coefficient <- function(x) {
    x <- sort(x)
    n <- length(x)
    G <- sum((2 * (1:n) - n - 1) * x) / (n * sum(x))
    return(1 - G)
  }
  
  # Calculating Simpson's Index
  simpsons_index <- function(x) {
    n <- sum(x)
    D <- sum((x / n) ^ 2)
    return(1 - D)
  }
  
  # Calculating Shannon Index
  shannon_index <- function(x) {
    p <- x / sum(x)
    H <- -sum(p * log(p))
    return(H)
  }
  
  # Access the nested list structure and unlist it
  refined_list <- unlist(unname(gene_lists))
  
  refined_list %>% unname() %>% unlist %>% table() %>% as.vector() -> frequency_genes
  
  # Calculate entropy, ensuring non-negative values if necessary
  # entropy_val <- if (all(refined_list >= 0)) entropy(refined_list) else NA
  
  # Calculate total genes
  total_genes <- length(refined_list)
  
  # Calculate unique gene count
  unique_genes <- length(unique(refined_list))
  
  gini <- gini_coefficient(frequency_genes)
  simpsons <- simpsons_index(frequency_genes)
  shannon <- shannon_index(frequency_genes)
  
  # Create a data frame with each metric as a column
  results_table <- data.frame(
    entropy = entropy(refined_list),
    total_genes = total_genes,
    unique_genes = unique_genes,
    gini_index = gini,
    simpsons_index = simpsons,
    shannon_index = shannon
  )
  
  return(results_table)
}

################################################################################
filter_vector <- function(vector, value, operation = "%in%") {
  # Check if the operation is valid
  valid_operations <- c("==", "!=", "%in%", "!%in%")
  if (!operation %in% valid_operations) {
    stop("Unsupported operation. Use '==', '!=', '%in%', or '!%in%'.")
  }
  
  # Apply the filtering based on the specified operation
  if (operation == "==") {
    filtered_vector <- vector[vector == value]
  } else if (operation == "!=") {
    filtered_vector <- vector[vector != value]
  } else if (operation == "%in%") {
    # Ensure `value` is a vector when using `%in%`
    if (!is.vector(value)) {
      stop("For '%in%' operation, the 'value' must be a vector.")
    }
    filtered_vector <- vector[vector %in% value]
  } else if (operation == "!%in%") {
    # Ensure `value` is a vector when using '!%in%'
    if (!is.vector(value)) {
      stop("For '!%in%' operation, the 'value' must be a vector.")
    }
    filtered_vector <- vector[!vector %in% value]
  }
  
  return(filtered_vector)
}

################################################################################
create_secondary_list <- function(data, 
                                  ...,
                                  gene_column = "hgnc_symbol",
                                  gene_to_remove = ""){
  data %>% 
    filter(...) %>%
    .[[gene_column]] %>% 
    convert_genes_to_vector() %>% 
    filter_vector(., value = gene_to_remove, "!%in%")
}

################################################################################
# Function to create master lists based on cumulative sum thresholds
create_master_gene_lists <- function(cumsum_thresholds, gene_distribution) {
  # Automatically generate names based on the cumulative sum thresholds
  names_list <- paste("master_cumsum", cumsum_thresholds, sep = "_")
  
  # Map function to apply create_secondary_list based on cumulative sum thresholds
  master_lists <- setNames(
    purrr::map(cumsum_thresholds, function(x) {
      create_secondary_list(gene_distribution, 
                            cumulative_sum > x)
    }),
    names_list
  )
  return(master_lists)
}


################################################################################
# Function to create rare gene lists based on frequency thresholds
create_rare_gene_lists <- function(freq_thresholds, gene_distribution) {
  # Automatically generate names based on the input thresholds
  names_list <- sapply(freq_thresholds, function(x) {
    paste("freq", paste(x, collapse = ""), sep = "_")
  })
  
  # Apply the create_secondary_list function to each frequency threshold
  rare_gene_list <- lapply(freq_thresholds, function(x) {
    create_secondary_list(gene_distribution, 
                          freq %in% x)
  })
  
  # Set names for each list element based on the generated names_list
  setNames(rare_gene_list, names_list)
}


################################################################################
filter_list_by_names <- function(lst, names_to_filter, operation = "%in%") {
  # Check if the operation is valid
  valid_operations <- c("==", "!=", "%in%", "!%in%")
  if (!operation %in% valid_operations) {
    stop("Unsupported operation. Use '==', '!=', '%in%', or '!%in%'.")
  }
  
  # Get the names of the list elements
  list_names <- names(lst)
  
  # Apply the filtering based on the specified operation
  if (operation == "==") {
    filtered_names <- list_names[list_names == names_to_filter]
  } else if (operation == "!=") {
    filtered_names <- list_names[list_names != names_to_filter]
  } else if (operation == "%in%") {
    filtered_names <- list_names[list_names %in% names_to_filter]
  } else if (operation == "!%in%") {
    filtered_names <- list_names[!list_names %in% names_to_filter]
  }
  
  # Extract elements from the list based on filtered names
  filtered_list <- lst[filtered_names]
  
  return(filtered_list)
}


  
