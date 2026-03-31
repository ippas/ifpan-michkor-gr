regulation_norm_master_lists$up_50_log2ratio_nListPlus2


jaccard_index_multiple_comparison(
  regulation_master_data = regulation_norm_master_lists,
  tissue_data = tissue_master_df,
  parameters_df = regulation_norm_master_df %>% filter(size_list == 50, rank_criterion == "fdr")
) -> jaccard_results2

jaccard_results2 %>% bind_rows(.id = "name") %>% column_to_rownames(var = "name") %>% pheatmap()

jaccard_index_per_tissue(gene_list =  regulation_norm_master_lists$up_50_log2ratio_nListPlus2$gene_list,
                         tissue_data = tissue_master_df,
                         regulation = "up",
                         size_list = 50,
                         rank_criterion = "log2ratio") 

regulation_master_withour_norm_lists %>% 
  bind_rows(.id = "label") %>% 
  separate(label, into = c("regulation", "metric", "size_list")) %>% 
  filter(metric == "log2ratio") %>% 
  filter(regulation == "up") %>% 
  .$hgnc_symbol %>%
  jaccard_index_per_tissue(gene_list =  .,
                         tissue_data = tissue_master_df,
                         regulation = "up",
                         size_list = 50,
                         rank_criterion = "log2ratio") 

regulation_master_lists %>% 
  bind_rows(.id = "label") %>% 
  separate(label, into = c("regulation", "metric", "size_list")) %>% 
  filter(metric == "log2ratio") %>% 
  filter(regulation == "down") %>% 
  .$hgnc_symbol %>%
  jaccard_index_per_tissue(gene_list =  .,
                           tissue_data = tissue_master_df,
                           regulation = "down",
                           size_list = 50,
                           rank_criterion = "log2ratio") 
