regulation_normalized_master_list <- list()


parameters_regulation_norm_master_df <- 
  data_frame(
  regulation = rep(c("up", "down"), 2) %>% rep(., 4) %>% rep(4),
  rank_criterion = rep(c("log2ratio", "fdr"), each = 2) %>% rep(., 4) %>% rep(4),
  size_list = rep(c(10, 25, 50, 100), each = 4) %>% rep(4),
  top_genes = rep(c(10, 25, 50, 100), each = 4) %>% rep(4),
  additional_adjustment = rep(c(0, 1, 2, 3), each = 16)
)

lapply(1:nrow(parameters_regulation_norm_master_df), function(i){
  print(parameters_regulation_norm_master_df[i, ])
  
  tryCatch({
    create_master_list_normalized_per_tissue(tissue_master_list = tissue_master_lists, 
                                             regulation = parameters_regulation_norm_master_df[i, "regulation"] %>% as.character(), 
                                             size_list = parameters_regulation_norm_master_df[i, "size_list"] %>% as.numeric(), 
                                             rank_criterion = parameters_regulation_norm_master_df[i, "rank_criterion"] %>% as.character(), 
                                             top_genes = parameters_regulation_norm_master_df[i, "top_genes"] %>% as.numeric(), 
                                             additional_adjustment = parameters_regulation_norm_master_df[i, "additional_adjustment"] %>% as.numeric()
    )
  }, error = function(e) NULL)
  
}) -> regulation_norm_master_lists

# set names for results
names(regulation_norm_master_lists) <- parameters_regulation_norm_master_df %>% 
  mutate(names = paste(regulation, size_list, rank_criterion, paste0("nListPlus", additional_adjustment), sep = "_")) %>% .$names 

regulation_norm_master_lists %>% 
  lapply(., function(x){x$master_df}) %>% 
  bind_rows(.id = "label") %>% 
  separate(label, into = c("regulation", "size_list", "rank_criterion", "normalization"), sep = "_") -> regulation_norm_master_df

regulation_norm_master_df %>% 
  write_tsv_xlsx(., tsv_file = "results/google-drive/master-lists/regulation-master-list-rank-score-normalize.tsv")

regulation_norm_master_df %>% 
  filter(size_list == 50, normalize == "n_list + 2") %>% 
  write_tsv_xlsx(., tsv_file = "results/google-drive/master-lists/regulation-master-list-rank-score-top50-norm-nListPlus2.tsv")



