# Define the function
create_master_genes_double_rank_score <- function(data_master_tissue, regulation, rank_criterion, size_list, n_top) {
  data_master_tissue %>%
    filter(regulation == {{regulation}}, 
           rank_criterion == {{rank_criterion}}, 
           size_list == {{size_list}}) %>%
    group_by(tissue) %>%
    nest() %>% 
    mutate(data = map(data, ~mutate(.x, rank_tissue = c(50:1)))) %>%
    unnest(cols = c(data)) %>%
    select(hgnc_symbol, rank_tissue) %>%
    group_by(hgnc_symbol) %>%
    nest() %>%
    mutate(sum_rank_tissue = map(data, ~sum(.x$rank_tissue))) %>%
    unnest(cols = c(sum_rank_tissue)) %>%
    arrange(desc(sum_rank_tissue)) %>%
    head(n_top) %>% 
    select(-data)
}


