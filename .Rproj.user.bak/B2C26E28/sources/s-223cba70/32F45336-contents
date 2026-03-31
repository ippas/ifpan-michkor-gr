
gene_list_log2ratio_top_50 %>% 
  select(-data) %>% 
  unnest(top_50_genes) %>% 
  ungroup() %>%
  group_by(simple_tissue, regulation) %>% 
  nest() %>% 
  mutate(n_list = map(data, ~ .x %>% .$label_regulation %>% unique %>% length)) %>% 
  mutate(n_genes_tissue = map(data, ~ .x %>% .$hgnc_symbol %>% unique %>% length)) %>% 
  unnest(c(n_list, n_genes_tissue)) %>% 
  unnest(data) %>% 
  group_by(regulation, simple_tissue, hgnc_symbol) %>% 
  nest() %>% 
  mutate(sum_rs = map(data, ~ .x$rank_score %>% sum)) %>% 
  unnest(sum_rs) %>% 
  group_by(regulation, simple_tissue) %>% 
  nest() %>% 
  mutate(top_50_rs = map(data, ~ .x %>% arrange(desc(sum_rs)) %>% head(50))) %>% 
  select(-data) %>% 
  unnest(top_50_rs) -> top50_genes_sum_rs


gene_list_log2ratio_top_50 %>% 
  select(-data) %>% 
  unnest(top_50_genes) %>%
  filter(simple_tissue %in% c("lung", "blood", "brain")) %>% 
  mutate(simple_cell = case_when(
    cell %in% c("ASM", "ASM-siCEBPO", "ASM-siCNTR") ~ "ASM",
    cell %in% c("A549", "HBE", "BEAS-2B", "pHBECs", 
                     "H1944", "H1975", "H2122") ~ "epithelial",
    cell %in% c("macrophages", "hMDM", "mBMDM", "THP-1", "Monocyte") ~ "macrophage-like",
    cell %in% c("Bcell", "NKcell", "Tcell", "NALM6", "REH-overexpression-GCR") ~ "lymphoid",
    
    # Brain cells (from experiment names)
    grepl("_mglia_", label_regulation) ~ "microglia",
    grepl("prefrontal-cortex|hippocampus", label_regulation) ~ "brain_neuronal",
    
    TRUE ~ NA_character_
  )) %>% 
  filter(!is.na(simple_cell)) %>% 
  mutate(simple_tissue = simple_cell) %>% 
  ungroup() %>%
  group_by(simple_tissue, regulation) %>% 
  nest() %>% 
  mutate(n_list = map(data, ~ .x %>% .$label_regulation %>% unique %>% length)) %>% 
  mutate(n_genes_tissue = map(data, ~ .x %>% .$hgnc_symbol %>% unique %>% length)) %>% 
  unnest(c(n_list, n_genes_tissue)) %>% 
  unnest(data) %>% 
  group_by(regulation, simple_tissue, hgnc_symbol) %>% 
  nest() %>% 
  mutate(sum_rs = map(data, ~ .x$rank_score %>% sum)) %>% 
  unnest(sum_rs) %>% 
  group_by(regulation, simple_tissue) %>% 
  nest() %>% 
  mutate(top_50_rs = map(data, ~ .x %>% arrange(desc(sum_rs)) %>% head(50))) %>% 
  select(-data) %>% 
  unnest(top_50_rs) -> cell_top50_genes_sum_rs
  

top50_genes_sum_rs

cell_top50_genes_sum_rs
  

