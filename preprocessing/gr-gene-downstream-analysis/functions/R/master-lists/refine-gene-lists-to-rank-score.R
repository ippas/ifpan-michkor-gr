filter_unique_gene_from_list <- function(data, genes_column, arrange_by){
  data %>% 
    group_by(label, !!!sym(genes_column)) %>% 
    nest %>%
    mutate(data2 = map(data, ~{.x %>% arrange(desc(arrange_by)) %>% head(1)})) %>% 
    select(-data) %>% unnest(data2)
}


refine_gene_lists_to_rank_score <- function(data, columns = NULL, keep_column = NULL, genes_column = "hgnc_symbol", size_threshold = 0, arrange_by) {
  
  Sys.time() -> start_time
  
  # creating gene lists
  data %>% 
    mutate(label = paste(!!!syms(columns), sep = "_")) %>%
    # select(c(!!!sym(genes_column), label)) %>% 
    filter_unique_gene_from_list(data = ., genes_column = genes_column, arrange_by = arrange_by) -> labeled_data
  
  labeled_data %>% 
    split(.$label, .[[genes_column]]) %>% 
    lapply(., function(x){x[[genes_column]]}) -> gene_lists
  
  
  # creating log2ratio lists
  labeled_data %>% 
    split(.$label, .[[arrange_by]]) %>% 
    lapply(., function(x){x[[arrange_by]]}) %>% 
    lapply(., as.numeric) -> tmp
  
  assign(arrange_by, tmp)
  
  print(  labeled_data  %>% head)

  # filtering and labeling the data
  data %>%
    mutate(label = paste(!!!syms(columns), sep = "_")) %>%
    select(label, !!!syms(columns), !!!(keep_column)) -> metadata
  
  # Counting genes in each list and preparing a count dataframe
  gene_counts <- sapply(gene_lists, length)
  count_df <- data.frame(label = names(gene_counts), gene_count = gene_counts, stringsAsFactors = FALSE)
  
  # Adding gene counts to metadata
  metadata <- metadata %>%
    distinct() %>%
    left_join(count_df, by = "label")
  
  # filter by size gene lists
  lists_to_keep <- metadata %>% 
    filter(gene_count >= size_threshold) %>% 
    .$label
  
  metadata <- metadata %>% filter(label %in% lists_to_keep)
  
  gene_lists <- gene_lists[lists_to_keep]
  
  list(metadata = metadata,
       gene_lists = gene_lists) -> results
  
  results[[arrange_by]] <- get(arrange_by)
  
  Sys.time() -> end_time
  
  print(end_time - start_time)
  
  return(results)
}
