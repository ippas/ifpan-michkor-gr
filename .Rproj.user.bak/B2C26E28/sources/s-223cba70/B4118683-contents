################################################################################
# jaccard index, for master lists up and down
jaccard_index(
  regulation_master_df %>% 
    filter(regulation == "up",
           size_list == 50,
           rank_criterion == "log2ratio") %>% 
    .$hgnc_symbol,
  regulation_master_df %>% 
    filter(regulation == "up",
           size_list == 50,
           rank_criterion == "fdr") %>% 
    .$hgnc_symbol)


jaccard_index(
  regulation_master_df %>% 
    filter(regulation == "down",
           size_list == 50,
           rank_criterion == "log2ratio") %>% 
    .$hgnc_symbol,
  regulation_master_df %>% 
    filter(regulation == "down",
           size_list == 50,
           rank_criterion == "fdr") %>% 
    .$hgnc_symbol)


################################################################################
# upregulated genes
tissue_master_df %>%
  .$tissue %>%
  unique() %>%
  {lapply(., function(x){
    jaccard_index(
      tissue_master_df %>% 
        filter(regulation == "up",
               size_list == 50,
               rank_criterion == "log2ratio",
               tissue == x) %>% 
        .$hgnc_symbol,
      tissue_master_df %>% 
        filter(regulation == "up",
               size_list == 50,
               rank_criterion == "fdr",
               tissue == x) %>% 
        .$hgnc_symbol)
  })} %>%
  setNames(tissue_master_df %>% .$tissue %>% unique()) %>% 
  bind_rows(.id = "tissue") %>% t  %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "tissue") -> jaccard_up_df

# downregulated genes
tissue_master_df %>%
  .$tissue %>%
  unique() %>%
  {lapply(., function(x){
    jaccard_index(
      tissue_master_df %>% 
        filter(regulation == "down",
               size_list == 50,
               rank_criterion == "log2ratio",
               tissue == x) %>% 
        .$hgnc_symbol,
      tissue_master_df %>% 
        filter(regulation == "down",
               size_list == 50,
               rank_criterion == "fdr",
               tissue == x) %>% 
        .$hgnc_symbol)
  })} %>%
  setNames(tissue_master_df %>% .$tissue %>% unique()) %>% 
  bind_rows(.id = "tissue") %>% t %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "tissue") -> jaccard_down_df

left_join(jaccard_up_df, jaccard_down_df, by  = "tissue") %>% 
  set_colnames(c("tissue", "upregulated", "downregulated")) %>% 
  write_tsv_xlsx(., tsv_file = "results/google-drive/master-tissue-lists/jaccard-tissues.tsv")



################################################################################
# Function to calculate Jaccard indices for a given regulation type and rank criterion
calculate_jaccard_indices <- function(master_data, regulation, criterion) {
  unique_tissues <- unique(tissue_master_df$tissue)
  
  jaccard_results <- lapply(unique_tissues, function(tissue2) {
    master_genes <- regulation_master_df %>%
      filter(regulation == regulation, size_list == 50, rank_criterion == criterion) %>%
      .$hgnc_symbol
    
    tissue_genes <- tissue_master_df %>%
      filter(regulation == regulation, size_list == 50, rank_criterion == criterion, tissue == tissue2) %>%
      .$hgnc_symbol
    
    jaccard_index(master_genes, tissue_genes)
  }) %>% 
    setNames(unique_tissues) %>% 
    bind_rows(.id = "tissue") %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "tissue")
  
  return(jaccard_results)
}

# Apply the function for each combination of regulation and rank criterion
results_log2ratio_up <- calculate_jaccard_indices(regulation_master_df, "up", "log2ratio")
results_log2ratio_down <- calculate_jaccard_indices(regulation_master_df, "down", "log2ratio")
results_fdr_up <- calculate_jaccard_indices(regulation_master_df, "up", "fdr")
results_fdr_down <- calculate_jaccard_indices(regulation_master_df, "down", "fdr")

# Combine all results into a single data frame for easy comparison or export
final_results <- list(
  log2ratio_up = results_log2ratio_up,
  log2ratio_down = results_log2ratio_down,
  fdr_up = results_fdr_up,
  fdr_down = results_fdr_down
)

final_results %>% bind_rows(.id = "comparison") %>% 
  rename(jaccard = V1) %>% 
  spread(., comparison, jaccard) %>% 
  column_to_rownames("tissue") %>% 
  as.matrix() %>% pheatmap()



################################################################################
# upregulated genes
tissue_master_df %>% filter(size_list == 50, regulation == "up", rank_criterion == "log2ratio") %>% .$tissue %>% unique() %>% 
  {lapply(., function(x){
    jaccard_index(
      tissue_master_df %>% 
        filter(regulation == "up",
               size_list == 50,
               rank_criterion == "log2ratio",
               tissue == x) %>% 
        .$hgnc_symbol,
      regulation_norm_master_lists$up_50_log2ratio_nListPlus2$gene_list)
  })} %>% 
  setNames(tissue_master_df %>% filter(size_list == 50, regulation == "up", rank_criterion == "log2ratio") %>%  .$tissue %>% unique()) %>% 
  bind_rows(.id = "tissue") %>% t  %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "tissue") %>% column_to_rownames(var = "tissue") %>% pheatmap()
