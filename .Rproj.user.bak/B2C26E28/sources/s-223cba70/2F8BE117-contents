################################################################################
# prepare data frame with parameters to prepare master lists for tissues
data.frame(tissue = papers_data_preprocessing %>% 
             filter(!(detailed_tissue %in% c("Other"))) %>% 
             .$detailed_tissue %>% unique() %>% rep(each = 8),
           metric = rep(c("log2ratio", "fdr"), 40),
           sort_order = rep(c("decrease", "increase"), 40),
           top = rep(rep(c(10, 25, 50, 100), each = 2), 10)
) -> parameters_df_detailed_tissues

rbind({parameters_df_detailed_tissues %>% mutate(regulation = "up")}, 
      {parameters_df_detailed_tissues %>% mutate(regulation = "down")}) %>% 
  mutate(sort_order = ifelse(regulation == "down", "increase", sort_order)) -> parameters_df_detailed_tissues

lapply(1:nrow(parameters_df_detailed_tissues), function(i){
  print(parameters_df_detailed_tissues[i, ])
  
  tryCatch({
    create_ranked_master_list(data = papers_data_preprocessing,
                              arrange_by = parameters_df_detailed_tissues[i, "metric"],
                              detailed_tissue == parameters_df_detailed_tissues[i, "tissue"],
                              regulation == parameters_df_detailed_tissues[i, "regulation"],
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
                              keep_column =  c("simple_tissue", "method", "detailed_tissue"),
                              size_list = parameters_df_detailed_tissues[i, "top"],
                              max_rank = parameters_df_detailed_tissues[i, "top"],
                              top_n = parameters_df_detailed_tissues[i, "top"],
                              sort_order = parameters_df_detailed_tissues[i, "sort_order"]
          
    )
  }, error = function(e) NULL)
  
}) -> detailed_tissue_master_lists

# set names for results
names(detailed_tissue_master_lists) <- parameters_df_detailed_tissues %>% mutate(names = paste(tissue, regulation, top, metric, sep = "_")) %>% .$names 


detailed_tissue_master_lists %>% lapply(., function(x){
  x$master_df
}) %>% bind_rows(.id = "label") %>% 
  separate(label, into = c("tissue", "regulation", "size_list", "rank_criterion"), sep = "_") -> detailed_tissue_master_df


compute_rank_score_log2ratioPlusFDR(data = detailed_tissue_master_lists,  n_top = 50) -> detailed_tissue_master_lists_rsLog2ratioPlusFDR

detailed_tissue_master_lists_rsLog2ratioPlusFDR %>% 
  bind_rows(.id = "name") %>% 
  separate(name, into = c("tissue", "regulation", "size_list", "rank_criterion"), sep = "_") %>% 
  rename(sum_score = score) %>% 
  rbind(., tissue_master_df) %>% group_by(tissue) %>% nest %>% 
  mutate(tissue_rank_score = map(data, ~{
    n <- nrow(.x)
    rank_scores <- seq(50, 50 - n + 1)
    .x <- .x %>% mutate(tissue_rank_score = rank_scores)
    .x
  })) %>%
  unnest(cols = c(data, tissue_rank_score)) -> detailed_tissue_master_df


################################################################################
# create master gene lists
create_master_genes_double_rank_score(
  data_master_tissue =  detailed_tissue_master_df,
  regulation = "up",
  rank_criterion = "log2ratio",
  n_top = 50,
  size_list = 50
) -> regulation_master_double_rank_score_lists$detailed_up_log2ratio_50

create_master_genes_double_rank_score(
  data_master_tissue =  detailed_tissue_master_df,
  regulation = "down",
  rank_criterion = "log2ratio",
  n_top = 50,
  size_list = 50
) -> regulation_master_double_rank_score_lists$detailed_down_log2ratio_50


################################################################################
# add general list for brain, blood, lung
detailed_tissue_master_df %>% 
  filter(!(tissue %in% c("brain", "blood", "lung"))) %>% 
  rbind(., tissue_master_df %>% 
          filter(tissue %in% c("brain", "blood", "lung"))
  ) %>% 
  mutate(tissue = case_when(
    tissue == "brain" ~ "brain-general",
    tissue == "lung" ~ "lung-general",
    tissue == "blood" ~ "blood-general", 
    TRUE ~ tissue
  )) -> detailed_tissue_master_df


################################################################################
# save gene lists to files
detailed_tissue_master_df %>% 
  write_tsv_xlsx(data = ., 
                 tsv_file= "results/google-drive/master-tissue-lists/detailed_tissue-master-list-rank-score.tsv")


detailed_tissue_master_df %>% 
  filter(size_list == 50) %>% 
  filter(rank_criterion == "log2ratio") %>% 
  write_tsv_xlsx(data = ., 
                 tsv_file= "results/google-drive/master-tissue-lists/detailed_tissue-master-list-rank-score-top50.tsv")

regulation_master_double_rank_score_lists[c("detailed_up_log2ratio_50", "detailed_down_log2ratio_50" )] %>% 
  bind_rows(.id = "label") %>% 
  separate(label, into = c("type", "regulation", "metric", "size_list")) %>%
  select(-type) %>% 
  write_tsv_xlsx(tsv_file = "results/google-drive/master-lists/reg-master-list-double-rank-score-normalzie-top50-detailed-tissue.tsv")


################################################################################
# overalp coeffient matrix
detailed_tissue_master_df %>%
  ungroup() %>%
  filter(size_list == 50, regulation == "up", rank_criterion == "log2ratio") %>%
  select(tissue, n_lists_filt) %>%
  unique() %>%
  mutate(name = paste0(tissue, " (", n_lists_filt, ")")) %>% select(-n_lists_filt) %>%  { set_names(.$name, .$tissue) }-> cf_heatmap_names


detailed_tissue_master_df %>% 
  filter(regulation == "up", 
         rank_criterion == "log2ratio",
         size_list == 50) %>% 
  ungroup(hgnc_symbol) %>%
  mutate(tissue = paste0(tissue, " (", n_lists_filt, ")")) %>% 
  select(c(tissue, hgnc_symbol)) %>% 
  as.data.frame() %>% 
  split(.$tissue, .$hgnc_symbol ) %>% lapply(., function(x) {x$hgnc_symbol}) %>%  
  append(list(master_up_double_rs = regulation_master_double_rank_score_lists$up_log2ratio_50$hgnc_symbol)) %>%
  overlap_coefficient_matrix(gene_list = .) -> cf_detailed_tissue_up_matrix 
  
svg("results/google-drive/figures/cf-heatmap-detailed-tissue-up.svg", width = 9, height = 9)
plot_triangular_heatmap(mat = cf_detailed_tissue_up_matrix, 
                        palette = as.vector(custom_palette), 
                        legend_name = "Overlap coefficient", 
                        scale_palette = TRUE)
dev.off()

  

detailed_tissue_master_df %>% 
  filter(regulation == "down", 
         rank_criterion == "log2ratio",
         size_list == 50) %>% 
  ungroup(hgnc_symbol) %>% 
  mutate(tissue = paste0(tissue, " (", n_lists_filt, ")")) %>% 
  select(c(tissue, hgnc_symbol)) %>%
  as.data.frame() %>% 
  split(.$tissue, .$hgnc_symbol ) %>% lapply(., function(x) {x$hgnc_symbol}) %>%  
  append(list(master_up_double_rs = regulation_master_double_rank_score_lists$down_log2ratio_50$hgnc_symbol)) %>%
  overlap_coefficient_matrix(gene_list = .) -> cf_detailed_tissue_down_matrix
  
dev.off()
svg("results/google-drive/figures/cf-heatmap-detailed-tissue-down.svg", width = 9, height = 9)
plot_triangular_heatmap(mat = cf_detailed_tissue_down_matrix, as.vector(custom_palette), legend_name = "Overlap coefficient", scale_palette = TRUE) 
dev.off()

wilcox.test(cf_detailed_tissue_up_matrix[lower.tri( cf_detailed_tissue_up_matrix)],
            cf_detailed_tissue_down_matrix[lower.tri( cf_detailed_tissue_down_matrix)])
