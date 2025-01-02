create_master_list_normalized_per_tissue <- function(tissue_master_list, regulation, size_list, rank_criterion, top_genes, additional_adjustment = 0){
  paste(regulation, size_list, rank_criterion, sep = "_") -> pattern
  
  tissue_master_lists[grepl(pattern, names(tissue_master_lists))] %>% 
    # Filter(function(x) !is.null(x), .)
    lapply(., function(x){
      x$genes_sum_score_df 
    }) %>% Filter(function(x) !is.null(x), .) %>% 
    lapply(., function(x){
      x %>% 
        mutate(n_lists_adj = n_lists_filt + additional_adjustment) %>% 
        mutate(tissue_norm_score = sum_score/n_lists_adj) 
    }) %>% 
    bind_rows(.id = "names") %>% 
    separate(names, into = c("tissue", "regulation", "size_list", "rank_criterion"), sep = "_") -> preprocessing_data
  
  summary <- list()
  
  print("hello")
  
  preprocessing_data %>% select(c(tissue, n_papers_raw))  %>%  unique() %>% .$n_papers_raw %>% summary() -> summary$n_papers_raw
  preprocessing_data %>% select(c(tissue, n_lists_raw))  %>%  unique() %>% .$n_lists_raw %>% summary() -> summary$n_lists_raw
  preprocessing_data %>% select(c(tissue, n_papers_filt))  %>%  unique() %>% .$n_papers_filt %>% summary() -> summary$n_papers_filt
  preprocessing_data %>% select(c(tissue, n_lists_filt))  %>%  unique() %>% .$n_lists_filt %>% summary() -> summary$n_lists_filt

  preprocessing_data %>% 
    group_by(hgnc_symbol) %>% 
    nest() %>% 
    mutate(sum_norm_score = map(data, ~sum(.x$tissue_norm_score))) %>%
    unnest(sum_norm_score) %>% 
    arrange(desc(sum_norm_score)) %>% 
    select(-data) %>% 
    head(as.numeric(top_genes)) %>% 
    mutate(normalize = paste0("n_list + ", additional_adjustment)) -> master_df
  
  print("hello")
  
  list(
    master_df = master_df,
    summary = summary,
    gene_list = master_df$hgnc_symbol
  ) -> results
  
  return(results)
  
}
