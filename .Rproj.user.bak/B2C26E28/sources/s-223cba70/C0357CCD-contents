top50_genes_sum_rs %>% 
  select(-data) %>% 
  ungroup %>% 
  select(hgnc_symbol, regulation, simple_tissue, sum_rs) %>% 
  write.table("results/figures/harmonized-tissue-gene-lists/summary-tissue-gene-lists/tables/gr-tissue-signatures.tsv", 
              sep = "\t",
              quote = F,
              row.names = F,
              col.names = T)


top50_genes_sum_rs %>% 
  select(-data) %>% 
  ungroup %>% 
  select(hgnc_symbol, regulation, simple_tissue, sum_rs) %>% 
  filter(simple_tissue %in% c("lung", "brain", "blood")) %>% 
  write.table("results/figures/harmonized-tissue-gene-lists/summary-tissue-gene-lists/tables/gr-tissue-signatures-lungBloodBrain.csv", 
              sep = ";",
              quote = F,
              row.names = F,
              col.names = T)
