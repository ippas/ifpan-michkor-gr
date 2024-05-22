refine_and_label_gene_lists <- function(data, columns = NULL, keep_column = NULL, genes_column = "hgnc_symbol", size_threshold = 10) {
  
  # creating gene lists
  data %>% 
    mutate(label = paste(!!!syms(columns), sep = "_")) %>%
    select(c(!!!sym(genes_column), label)) %>% 
    split(.$label, .[[genes_column]]) %>% 
    # lapply(., unique) %>%
    lapply(., function(x){x[[genes_column]]}) -> gene_lists
  
  
  # creating log2ratio lists
  data %>% 
    mutate(label = paste(!!!syms(columns), sep = "_")) %>%
    select(c(log2ratio, label)) %>% 
    split(.$label, .[["log2ratio"]]) %>% 
    # lapply(., unique) %>% 
    lapply(., function(x){x[["log2ratio"]]}) %>% 
    lapply(., as.numeric) -> log2ratio_lists
  
  
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
       gene_lists = gene_lists, 
       log2ratio_lists = log2ratio_lists) -> data
  
  return(data)
}
