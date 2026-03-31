compute_ranked_scores_for_genes <- function(data, max_rank = 50, rank_criterion = "log2ratio", sort_order = "decrease"){
  metadata <- data$metadata
  
  calculate_ranked_scores <- function(data, gene_list_name, rank_criterion, sort_order){
    
    gene_lists <- data$gene_lists
    data_to_calculate_rank <- data[[rank_criterion]]
    
    # Create a data frame
    ranked_df <- data.frame(
      hgnc_symbol = gene_lists[[gene_list_name]],
      value = data_to_calculate_rank[[gene_list_name]]
    )
    
    # Sort based on the direction specified by the sort_order argument
    if (sort_order == "decrease") {
      ranked_df <- ranked_df %>%
        dplyr::arrange(desc(value))  # Sort in descending order
    } else if (sort_order == "increase") {
      ranked_df <- ranked_df %>%
        dplyr::arrange(value)  # Sort in ascending order
    }
    
    # Calculate rank and score
    ranked_df <- ranked_df %>%
      dplyr::mutate(
        rank = dplyr::row_number(),  # Adds a new column 'rank' with values from 1 to n
        score = pmax(max_rank + 1 - rank, 0)  # Calculates the score
      )
    
    # Set column names
    names(ranked_df) <- c("hgnc_symbol", rank_criterion, "rank", "score")
    
    return(ranked_df)
  }
  
  print(metadata)
  ranked_scores_data <- lapply(metadata$label, function(x){
    calculate_ranked_scores(data = data, gene_list_name = x, rank_criterion = rank_criterion, sort_order = sort_order)
  })
  
  names(ranked_scores_data) <- metadata$label
  
  data$ranked_scores_data <- ranked_scores_data
  
  return(data)
}


# Define the function
calculate_top_genes <- function(data, top_n = 50) {
  data %>%
    bind_rows(., .id = "list_name") %>%  # Bind rows and add an identifier column
    group_by(hgnc_symbol) %>%           # Group by gene symbol
    nest() %>%                          # Nest the data within each group
    mutate(sum_score = map(data, ~ sum(.x$score))) %>%  # Calculate sum of scores
    unnest(sum_score) %>%               # Unnest the sum_score column
    arrange(desc(sum_score)) %>%        # Arrange in descending order of sum_score
    select(hgnc_symbol, sum_score) %>%  # Select the necessary columns
    head(top_n) %>%                     # Get the top 'n' entries
    pull(hgnc_symbol)                   # Extract only the hgnc_symbol column
}


################################################################################
create_ranked_master_list <-
  function(data,
           ...,
           # treatment_to_remove,
           arrange_by,
           top_n = 50,
           max_rank = 50,
           size_list = 0,
           columns,
           keep_column,
           sort_order = "decrease") {
    
    args <- enquos(...)
    print(args) 
    
    data %>%
      filter(!!!args) -> data
    
    print(dim(data))
    
    data %>%
      drop_na(!!!sym(arrange_by)) %>%
      filter(get(arrange_by) != "NA") %>%
      refine_gene_lists_to_rank_score(
        data = .,
        columns = columns,
        size_threshold = 0,
        # keep_column =  c("simple_tissue", "method"),
        keep_column = {{keep_column}},
        arrange_by = arrange_by,
      ) -> preprocessing_data
    
    print(preprocessing_data)
    
    preprocessing_data$metadata$source %>% unique() %>% length() -> n_papers_raw
    preprocessing_data$metadata$label %>% unique()  %>% length() -> n_lists_raw
    preprocessing_data %>%  filter_gene_list(gene_count > size_list, log2ratio = TRUE) -> preprocessing_data
    
    preprocessing_data$metadata$source %>% unique() %>% length() -> n_papers_filtered
    preprocessing_data$metadata$label %>% unique()  %>% length() -> n_lists_filtered
    
    compute_ranked_scores_for_genes(data = preprocessing_data,
                                    max_rank = max_rank,
                                    rank_criterion = arrange_by,
                                    sort_order = sort_order) -> ranked_data
    
    print(ranked_data)
    
    ranked_data$ranked_scores_data %>%
      lapply(., function(x){x %>% filter(score > 0)}) %>%
      lapply(., function(x){x$hgnc_symbol}) -> gene_lists_top

    
    ranked_data$ranked_scores_data %>%
      bind_rows(., .id = "list_name") %>%
      group_by(hgnc_symbol) %>%
      nest() %>%
      mutate(sum_score = map(data, ~sum(.x$score))) %>%
      unnest(sum_score) %>%
      arrange(desc(sum_score)) %>%
      select(c(hgnc_symbol, sum_score)) %>% 
      mutate(n_papers_raw = n_papers_raw,
             n_lists_raw = n_lists_raw,
             n_papers_filt = n_papers_filtered,
             n_lists_filt = n_lists_filtered
      ) -> genes_sum_score_df
    
 
    
    genes_sum_score_df %>% 
      head(top_n) %>%
      mutate(n_papers_raw = n_papers_raw,
             n_lists_raw = n_lists_raw,
             n_papers_filt = n_papers_filtered,
             n_lists_filt = n_lists_filtered
      ) -> master_df
    

    
    ranked_data$ranked_metric <- arrange_by
    ranked_data$gene_lists_top <- gene_lists_top
    ranked_data$genes_sum_score_df <- genes_sum_score_df
    ranked_data$master_df <- master_df
    
    print("#################################################################")
    print(ranked_data)
    
    return(ranked_data)
    
  }

