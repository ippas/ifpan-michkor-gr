papers_data_preprocessing %>% 
  filter(simple_tissue == "brain") %>% 
  mutate(label2 = paste(label, regulation, sep = "_")) %>% 
  filter(regulation == "down") %>% 
  .$label2 %>%  unique %>% length()

papers_data_preprocessing %>% 
  filter(simple_tissue == "brain") %>% 
  mutate(label2 = paste(label, regulation, sep = "_")) %>% 
  .$label2 %>% 
  table %>% as.data.frame() %>% .$Freq %>% summary


papers_data_preprocessing %>% 
  filter(simple_tissue == "brain") %>% 
  mutate(label2 = paste(label, regulation, sep = "_")) %>% 
  .$hgnc_symbol %>% table %>% 
  as.data.frame() %>% 
  .$Freq %>% table

library(e1071)
skewness(
  papers_data_preprocessing %>% 
    filter(simple_tissue == "brain") %>% 
    mutate(label2 = paste(label, regulation, sep = "_")) %>% 
    .$hgnc_symbol %>% 
    table %>% 
    as.data.frame() %>% 
    .$Freq
)

papers_data_preprocessing %>% 
  filter(simple_tissue == "brain") %>% 
  mutate(label2 = paste(label, regulation, sep = "_")) %>% 
  group_by(label2) %>% 
  nest() %>% 
  mutate(n_genes = map(data, ~ .x %>% nrow)) %>% 
  unnest() %>% 
  # filter(n_genes > 10) %>% 
  ungroup() %>% 
  group_by(hgnc_symbol) %>% 
  nest() %>% 
  mutate(freq_gene = map(data, ~ .x %>% nrow)) %>% 
  unnest() %>% 
  ungroup %>%
  as.data.frame() %>% 
  filter(freq_gene == 10) %>%
  # dplyr::select(regulation, hgnc_symbol) %>% 
  unique %>% 
  group_by(hgnc_symbol) %>% 
  nest() %>% 
  mutate(n_reg = map(data, ~ .x %>% nrow)) %>% 
  unnest() %>% 
  # filter(regulation == "up") %>% 
  filter(n_reg == 2) %>% 
  .$hgnc_symbol


papers_data_preprocessing %>% 
  filter(simple_tissue == "brain") %>% 
  mutate(label2 = paste(label, regulation, sep = "_")) %>% 
  filter(is.na(log2ratio)) %>% 
  .$label2 %>% 
  table %>% as.data.frame() %>% 
  filter(Freq >= 50)

papers_data_preprocessing %>% 
  filter(simple_tissue == "brain") %>% head

papers_data_preprocessing %>% 
  filter(simple_tissue == "brain")



metaphenotypes_AllBrain2Brain2Global$overlap$original_data$grouped_data %>%
  filter(Var2 %in% c("alcohol_use_and_misuse", "anxiety_and_nervousness",
                     "cognition_and_processing_speed", "depressive_symptomatology",
                     "trauma")) %>%
  select(-c(n_unique_genes, n_unique_genes_p0.05, n_phenotypes_genes)) %>%
  unnest() %>%
  filter(p_value < 0.05) %>%
  as.data.frame()
  pull(overlap_genes) %>%
  strsplit(",") %>%
  unlist() %>%
  cat(sep = "\n")
  
  
papers_data_preprocessing %>% 
  filter(simple_tissue == "brain") %>% 
  mutate(label2 = paste(label, regulation, sep = "_")) %>% 
  filter(!is.na(log2ratio)) %>% 
  group_by(label2) %>% 
  nest() %>% 
  mutate(n_geneInList = map(data, ~ nrow(.x))) %>% 
  unnest(n_geneInList) %>% 
  filter(n_geneInList >= 50) %>% 
  unnest(data) %>% 
  as.data.frame() %>% 
  group_by(label2)



papers_data_preprocessing %>% 
  filter(simple_tissue == "brain") %>% 
  mutate(label2 = paste(label, regulation, sep = "_"),
         log2ratio = as.numeric(log2ratio)) %>% 
  filter(!is.na(log2ratio)) %>% 
  group_by(label2) %>% 
  nest() %>% 
  mutate(n_geneInList = map_int(data, nrow)) %>% 
  filter(n_geneInList >= 50) %>% 
  unnest(data) %>% 
  group_by(label2) %>% 
  arrange(desc(abs(log2ratio)), .by_group = TRUE) %>% 
  slice_head(n = 50) %>% 
  mutate(rank = 51 - row_number()) %>%   # największe dostaje 50, najmniejsze 1
  ungroup() %>% 
  as.data.frame() %>% 
  select(label2, hgnc_symbol, source, log2ratio, rank) %>% 
  filter(hgnc_symbol %in% c("FMO2", "CYP1B1", "PRG4", "FLNC", "IL1R1"))
