# prepare data frame with parameters to prepare master lists for tissues
data.frame(tissue = papers_data_preprocessing %>% 
             filter(!(simple_tissue %in% c("Other"))) %>% 
             .$simple_tissue %>% unique() %>% rep(each = 2),
           metric = rep(c("log2ratio", "fdr"), 14)
) %>% mutate(sort_order = ifelse(metric == "log2ratio", "decrease", "increase")) -> parameters_df_tissues_noReg



lapply(1:nrow(parameters_df_tissues_noReg), function(i){
  print(parameters_df_tissues_noReg[i, ])
  
  tryCatch({
    create_ranked_master_list(data = papers_data_preprocessing %>% mutate(log2ratio = abs(as.numeric(log2ratio))),
                              arrange_by = parameters_df_tissues_noReg[i, "metric"],
                              simple_tissue == parameters_df_tissues_noReg[i, "tissue"],
                              regulation %in% c("up", "down"),
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
                              size_list = 50,
                              max_rank = 50,
                              top_n = 100,
                              sort_order = parameters_df_tissues_noReg[i, "sort_order"]
    )
  }, error = function(e) NULL)
  
}) -> tissue_master_lists_noReg

# set names for results
names(tissue_master_lists_noReg) <- parameters_df_tissues_noReg %>% mutate(names = paste(tissue, "100", metric, sep = "_")) %>% .$names 


tissue_master_lists_noReg %>% lapply(., function(x){
  x$master_df 
}) %>% bind_rows(.id = "label") %>% 
  separate(label, into = c("tissue", "size_list", "rank_criterion"), sep = "_") -> tissue_master_noReg_df

write_tsv_xlsx(tissue_master_noReg_df, tsv_file = "results/google-drive/master-tissue-lists/tissue-master-100-noReg.tsv")



