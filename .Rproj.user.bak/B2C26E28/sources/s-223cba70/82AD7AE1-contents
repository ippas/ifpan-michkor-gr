################################################################################
# prepare data frame with parameters to prepare master lists for tissues
data.frame(tissue = papers_data_preprocessing %>% 
             filter(!(simple_tissue %in% c("Other"))) %>% 
             .$simple_tissue %>% unique() %>% rep(each = 8),
           metric = rep(c("log2ratio", "fdr"), 56),
           sort_order = rep(c("decrease", "increase"), 56),
           top = rep(rep(c(10, 25, 50, 100), each = 2), 14)
           ) -> parameters_df_tissues

rbind({parameters_df_tissues %>% mutate(regulation = "up")}, 
      {parameters_df_tissues %>% mutate(regulation = "down")}) %>% 
  mutate(sort_order = ifelse(regulation == "down", "increase", sort_order)) -> parameters_df_tissues

lapply(1:nrow(parameters_df_tissues), function(i){
  print(parameters_df_tissues[i, ])
  
  tryCatch({
    create_ranked_master_list(data = papers_data_preprocessing,
                              arrange_by = parameters_df_tissues[i, "metric"],
                              simple_tissue == parameters_df_tissues[i, "tissue"],
                              regulation == parameters_df_tissues[i, "regulation"],
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
                              size_list = parameters_df_tissues[i, "top"],
                              max_rank = parameters_df_tissues[i, "top"],
                              top_n = parameters_df_tissues[i, "top"],
                              sort_order = parameters_df_tissues[i, "sort_order"]
    )
  }, error = function(e) NULL)
  
}) -> tissue_master_lists

# set names for results
names(tissue_master_lists) <- parameters_df_tissues %>% mutate(names = paste(tissue, regulation, top, metric, sep = "_")) %>% .$names 



tissue_master_lists %>% lapply(., function(x){
  x$master_df
}) %>% bind_rows(.id = "label") %>% 
  separate(label, into = c("tissue", "regulation", "size_list", "rank_criterion"), sep = "_") -> tissue_master_df


compute_rank_score_log2ratioPlusFDR(data = tissue_master_lists,  n_top = 50) -> tissue_master_lists_rsLog2ratioPlusFDR

tissue_master_lists_rsLog2ratioPlusFDR %>% 
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
  unnest(cols = c(data, tissue_rank_score)) -> tissue_master_df


tissue_master_df %>% 
  write_tsv_xlsx(data = ., 
               tsv_file= "results/google-drive/master-tissue-lists/tissue-master-list-rank-score.tsv")


tissue_master_df %>% 
  filter(size_list == 50) %>% 
  write_tsv_xlsx(data = ., 
                 tsv_file= "results/google-drive/master-tissue-lists/tissue-master-list-rank-score-top50.tsv")


tissue_master_df %>% 
  filter(size_list == 50) %>% 
  filter(rank_criterion == "log2ratio") %>% 
  write_tsv_xlsx(data = ., 
                 tsv_file= "results/google-drive/master-tissue-lists/tissue-master-list-rs-log2ratio-top50.tsv")

################################################################################
regulation_master_double_rank_score_lists <- list()

create_master_genes_double_rank_score(data_master_tissue = tissue_master_df, regulation = "up", rank_criterion = "log2ratio", n_top = 50, size_list = 50) -> regulation_master_double_rank_score_lists$up_log2ratio_50
create_master_genes_double_rank_score(data_master_tissue = tissue_master_df, regulation = "up", rank_criterion = "fdr", n_top = 50, size_list = 50) -> regulation_master_double_rank_score_lists$up_fdr_50
create_master_genes_double_rank_score(data_master_tissue = tissue_master_df, regulation = "down", rank_criterion = "log2ratio", n_top = 50, size_list = 50) -> regulation_master_double_rank_score_lists$down_log2ratio_50
create_master_genes_double_rank_score(data_master_tissue = tissue_master_df, regulation = "down", rank_criterion = "fdr", n_top = 50, size_list = 50) -> regulation_master_double_rank_score_lists$down_fdr_50
create_master_genes_double_rank_score(data_master_tissue = tissue_master_df, regulation = "up", rank_criterion = "log2ratioPlusFDR", n_top = 50, size_list = 50)-> regulation_master_double_rank_score_lists$up_log2ratioPlusFDR_50
create_master_genes_double_rank_score(data_master_tissue = tissue_master_df, regulation = "down", rank_criterion = "log2ratioPlusFDR", n_top = 50, size_list = 50) -> regulation_master_double_rank_score_lists$down_log2ratioPlusFDR_50


regulation_master_double_rank_score_lists %>% 
  bind_rows(.id = "label") %>% 
  separate(label, into = c("regulation", "metric", "size_list")) %>% 
  write_tsv_xlsx(tsv_file = "results/google-drive/master-lists/reg-master-list-double-rank-score-normalzie-top50.tsv")

regulation_master_double_rank_score_lists %>% 
  bind_rows(.id = "label") %>% 
  separate(label, into = c("regulation", "metric", "size_list")) %>% 
  filter(metric == "log2ratio") %>% 
  write_tsv_xlsx(tsv_file = "results/google-drive/master-lists/reg-master-lists-double-rs-norm-log2ratio-top50.tsv")

################################################################################
results_enrichment_master_genes <- list()

regulation_master_lists %>% 
  bind_rows(.id = "label") %>% 
  separate(label, into = c("regulation", "metric", "size_list")) %>% 
  filter(metric == "log2ratio") %>% 
  filter(regulation == "up") %>% 
  .$hgnc_symbol %>%
  perform_enrichment_analysis(gene_list = , databases, fdr_threshold, overlap_threshold) -> results_enrichment_master_genes$master_up50_norm

regulation_master_lists %>% 
  bind_rows(.id = "label") %>% 
  separate(label, into = c("regulation", "metric", "size_list")) %>% 
  filter(metric == "log2ratio") %>% 
  filter(regulation == "down") %>% 
  .$hgnc_symbol %>% 
  perform_enrichment_analysis(gene_list = , databases, fdr_threshold, overlap_threshold) -> results_enrichment_master_genes$master_down50_norm

results_enrichment_master_genes$master_up50_norm %>% 
filter(enrichr_database == "MGI_Mammalian_Phenotype_Level_4_2021") %>%
  plotEnrich(., showTerms = 40, numChar = 50, y = "Count", orderBy = "Adjusted.P.value", title = "Upregulated")

results_enrichment_master_genes$master_down50_norm %>% 
  filter(enrichr_database == "MGI_Mammalian_Phenotype_Level_4_2021") %>% 
  plotEnrich(., showTerms = 40, numChar = 50, y = "Count", orderBy = "Adjusted.P.value", title = "Downregulated")


results_enrichment_master_genes$master_up50_norm %>% 
  filter(enrichr_database == "MGI_Mammalian_Phenotype_Level_4_2021") %>% nrow

results_enrichment_master_genes$master_up50_norm %>% 
  filter(enrichr_database == "MGI_Mammalian_Phenotype_Level_4_2021") %>% 
  .$Genes %>% convert_genes_to_vector(split = ";") %>% unique() %>% length()


results_enrichment_master_genes$master_down50_norm %>% 
  filter(enrichr_database == "MGI_Mammalian_Phenotype_Level_4_2021") %>% nrow


results_enrichment_master_genes$master_down50_norm %>% 
  filter(enrichr_database == "MGI_Mammalian_Phenotype_Level_4_2021") %>%
  .$Genes %>% convert_genes_to_vector(split = ";") %>% unique() %>% length()
