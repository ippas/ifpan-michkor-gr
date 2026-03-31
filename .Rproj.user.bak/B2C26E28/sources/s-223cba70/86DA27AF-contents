################################################################################
# prepare data frame with parameters to prepare master lists for 
parameters_df_tissues %>% 
  select(-tissue) %>% 
  unique() -> parameters_df_regulation


lapply(1:nrow(parameters_df_regulation), function(i){
  print(parameters_df_regulation[i, ])
  
  tryCatch({
    create_ranked_master_list(data = papers_data_preprocessing,
                              arrange_by = parameters_df_regulation[i, "metric"],
                              regulation == parameters_df_regulation[i, "regulation"],
                              !(treatment %in% treatment_to_remove),
                              columns = c(
                                "source",
                                "tissue",
                                "cell",
                                "dose",
                                "treatment",
                                "treatment_type",
                                "regulation",
                                "comparison",
                                "environment"
                              ),
                              keep_column =  c("simple_tissue", "method"),
                              size_list = parameters_df_regulation[i, "top"],
                              max_rank = parameters_df_regulation[i, "top"],
                              top_n = parameters_df_regulation[i, "top"],
                              sort_order = parameters_df_regulation[i, "sort_order"]
    )
  }, error = function(e) NULL)
  
}) -> regulation_master_without_norm_lists
# set names for results
names(regulation_master_without_norm_lists) <- parameters_df_regulation %>% mutate(names = paste(regulation, top, metric, sep = "_")) %>% .$names 


regulation_master_without_norm_lists %>% 
  lapply(., function(x){x$master_df}) %>% 
  bind_rows(.id = "name") %>% 
  separate(name, into = c("regulation", "size_list", "rank_criterion")) -> regulation_master_without_norm_df
