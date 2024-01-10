papers_data_preprocessing %>% 
  filter(!(source %in% c("pmid:NA", "marpiech_tissues_dex"))) %>% 
  .$source %>% unique()


papers_data_preprocessing %>% 
  filter(!(source %in% c("pmid:NA", "marpiech_tissues_dex"))) %>% 
  .$label %>% unique %>% length()

papers_data_preprocessing %>% 
  filter(!(source %in% c("pmid:NA", "marpiech_tissues_dex"))) %>% 
  select(c(source, tissue, cell)) %>% unique() %>% 
  write.table("results/presentation-summary-2023/pmid-tisue-cell.tsv", sep = "\t", 
              col.names = T, row.names = F)
  
terms_vector <- c("cortex", "embryos", "adrenal-gland", "hypothalamus", "lung", 
                  "embryos", "hippocampus", "lung", "brain", "brain", 
                  "brain", "placenta", "cortex", "hippocampus", "striatum", 
                  "striatum", "brain", "hippocampus", "hypothalamus", 
                  "hippocampus", "hippocampus", "hippocampus", "hippocampus", 
                  "cortex", "hippocampus", "hippocampus", "brain", "blood", 
                  "adipocyte-tissue", "adipocyte-tissue", "skeletal-muscle", 
                  "blood", "knee-cartilage", "skeletal-muscle", "lung", 
                  "liver", "blood", "kidney", "brain", "embryos", "kidney")

terms_vector %>% table %>% as.data.frame() %>% set_colnames(c("tissue", "number")) %>% 
  arrange(desc(number)) %>% 
  write.table("results/presentation-summary-2023/summary-tissue.tsv", sep = "\t", col.names = T, row.names = F)
