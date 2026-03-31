papers_data_preprocessing %>% 
  drop_na(log2ratio) %>%
  filter(log2ratio != "NA") %>% 
  filter(!(treatment %in% treatment_to_remove)) %>%
  filter(regulation == "up") %>%
  refine_gene_lists(data = .,
                    columns = c("source", "tissue", "cell", "dose", "treatment", "treatment_type", "regulation", "comparison", "environment"),
                    cumsum_thresholds = cumsum_thresholds,
                    freq_thresholds = freq_thresholds,
                    keep_column =  c("simple_tissue", "method")) -> tmp 


tmp$refine_gene_lists %>%  
  filter_gene_list(gene_count > 50, log2ratio = TRUE) -> tmp

compute_ranked_scores_for_genes(data = tmp) -> tmp

tmp$ranked_scores_data %>% 
  bind_rows(., .id = "list_name") %>% 
  group_by(hgnc_symbol) %>% 
  nest() %>% 
  mutate(sum_score = map(data, ~sum(.x$score))) %>% 
  unnest(sum_score) %>% 
  arrange(desc(sum_score)) %>% 
  select(c(hgnc_symbol, sum_score)) %>% head(50) %>% 
  rename(sum_ranked_score = sum_score) %>% 
  mutate(regulation = "up") -> gr_master_up


papers_data_preprocessing %>% 
  drop_na(log2ratio) %>%
  filter(log2ratio != "NA") %>% 
  filter(!(treatment %in% treatment_to_remove)) %>%
  filter(regulation == "down") %>%
  refine_gene_lists(data = .,
                    columns = c("source", "tissue", "cell", "dose", "treatment", "treatment_type", "regulation", "comparison", "environment"),
                    cumsum_thresholds = cumsum_thresholds,
                    freq_thresholds = freq_thresholds,
                    keep_column =  c("simple_tissue", "method")) -> tmp 


tmp$refine_gene_lists %>%  
  filter_gene_list(gene_count > 50, log2ratio = TRUE) -> tmp

compute_ranked_scores_for_genes(data = tmp) -> tmp

tmp$ranked_scores_data %>% 
  bind_rows(., .id = "list_name") %>% 
  group_by(hgnc_symbol) %>% 
  nest() %>% 
  mutate(sum_score = map(data, ~sum(.x$score))) %>% 
  unnest(sum_score) %>% 
  arrange(desc(sum_score)) %>% 
  select(c(hgnc_symbol, sum_score)) %>% head(50) %>% 
  rename(sum_ranked_score = sum_score) %>% 
  mutate(regulation = "down") -> gr_master_down




rbind(gr_master_up, gr_master_down) %>% 
  write_tsv_xlsx(tsv_file = "results/google-drive/master-lists/up-down-top50-mrs-size50.tsv")

rbind(gr_master_up, gr_master_down)
  write.table(file = "results/google-drive/master-lists//master-lists-up-down-ranked-score-top50.tsv", 
              sep = "\t", 
              row.names = F,
              col.names = T,
              quote = F)
  
write_xlsx(data, "results/summary-table-lists/master-lists-up-down-ranked-score-top50.xlsx")



