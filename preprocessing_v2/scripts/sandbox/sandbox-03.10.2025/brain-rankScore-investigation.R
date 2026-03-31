
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
  select(label2, hgnc_symbol, source, log2ratio, regulation, rank) %>% 
  # filter(regulation == "down") %>% 
  # filter(hgnc_symbol %in% c("FMO2", "CYP1B1", "PRG4", "FLNC", "IL1R1")) %>%
  group_by(hgnc_symbol, regulation) %>% 
  summarise(sum_rank = sum(rank), .groups = "drop") %>% 
  ungroup %>% 
  arrange(desc(sum_rank)) %>% 
  head(50) 


papers_data_preprocessing %>% 
  filter(simple_tissue == "brain") %>% 
  mutate(label2 = paste(label, regulation, sep = "_"),
         log2ratio = as.numeric(log2ratio)) %>% 
  filter(!is.na(log2ratio)) %>% 
  select(-index) %>% 
  select(-info) %>% 
  unique %>% 
  group_by(label2) %>% 
  nest() %>% 
  mutate(n_geneInList = map_int(data, nrow)) %>% 
  filter(n_geneInList >= 50) %>% 
  # w każdej liście robimy top50 + rank
  mutate(
    top50 = map(data, ~ .x %>% 
                  arrange(desc(abs(log2ratio))) %>% 
                  slice_head(n = 50) %>% 
                  mutate(rank = 51 - row_number()))
  ) %>% 
  select(label2, top50) %>% 
  unnest(top50) %>%
  # teraz globalne podsumowanie
  group_by(hgnc_symbol, regulation) %>% 
  summarise(
    sum_rank = sum(rank),
    n_lists  = n()   # w ilu listach gen się pojawił w top50
    # .groups = "drop"
  ) %>% 
  # filter(hgnc_symbol %in% c("CYP1B1", "FLNC", "FMO2", "IL1R1", "PRG4"))
  arrange(desc(sum_rank))
  filter(regulation == "up") %>%
  head(50) %>%   filter(hgnc_symbol %in% c("CYP1B1", "FLNC", "FMO2", "IL1R1", "PRG4"))


papers_data_preprocessing %>% 
  filter(simple_tissue == "brain") %>% 
  mutate(label2 = paste(label, regulation, sep = "_"),
         log2ratio = as.numeric(log2ratio)) %>% 
  filter(!is.na(log2ratio)) %>% 
  select(-index, -info) %>% 
  unique() %>% 
  group_by(label2) %>% 
  nest() %>% 
  mutate(n_geneInList = map_int(data, nrow)) %>% 
  filter(n_geneInList >= 50) %>% 
  # w każdej liście robimy top50 + rank
  mutate(
    top50 = map(data, ~ .x %>% 
                  arrange(desc(abs(log2ratio))) %>% 
                  slice_head(n = 50) %>% 
                  mutate(rank = 51 - row_number()))
  ) %>% 
  select(label2, top50) %>% 
  unnest(top50) %>%
  # globalne podsumowanie
  group_by(hgnc_symbol, regulation) %>% 
  summarise(
    sum_rank = sum(rank),
    n_lists  = n(),   
    lists    = list(cur_data()),   # tu nestujemy pełne dane źródłowe
    .groups = "drop"
  ) %>% 
  filter(regulation == "down") %>%
  arrange(desc(sum_rank)) %>% 
  head(50) %>% 
  filter(hgnc_symbol %in% c("CYP1B1", "FLNC", "FMO2", "IL1R1", "PRG4")) %>% unnest()
  
papers_data_preprocessing %>% 
  filter(simple_tissue == "brain") %>% 
  mutate(label2 = paste(label, regulation, sep = "_"),
         log2ratio = as.numeric(log2ratio)) %>% 
  filter(!is.na(log2ratio)) %>% 
  select(-index, -info) %>% 
  unique() %>% 
  group_by(label2) %>% 
  nest() %>% 
  mutate(n_geneInList = map_int(data, nrow)) %>% 
  filter(n_geneInList >= 50) %>% 
  mutate(
    top50 = map(data, ~ .x %>% 
                  arrange(desc(abs(log2ratio))) %>% 
                  slice_head(n = 50) %>% 
                  mutate(rank = 51 - row_number()))
  ) %>% 
  select(label2, top50) %>% 
  unnest(top50) %>%
  group_by(hgnc_symbol, regulation) %>% 
  summarise(
    sum_rank = sum(rank),
    n_lists  = n(),
    lists    = list(
      tibble(label2 = label2, log2ratio = log2ratio, rank = rank, source = source)
    ),
    .groups = "drop"
  ) %>%   filter(hgnc_symbol %in% c("CYP1B1", "FLNC", "FMO2", "IL1R1", "PRG4")) %>% 
  filter(regulation == "down") %>% unnest


papers_data_preprocessing %>% 
  filter(simple_tissue == "brain") %>% 
  mutate(label2 = paste(label, regulation, sep = "_"),
         log2ratio = as.numeric(log2ratio)) %>% 
  filter(!is.na(log2ratio)) %>% 
  select(-index, -info) %>% 
  unique() %>% 
  group_by(label2) %>% 
  nest() %>% 
  mutate(n_geneInList = map_int(data, nrow)) %>% 
  filter(n_geneInList >= 50) %>% 
  mutate(
    top50 = map(data, ~ .x %>% 
                  arrange(desc(abs(log2ratio))) %>% 
                  slice_head(n = 50) %>% 
                  mutate(rank = 51 - row_number()))
  ) %>% 
  select(label2, top50) %>% 
  unnest(top50) %>%
  group_by(hgnc_symbol, regulation) %>% 
  summarise(
    sum_rank = sum(rank),
    n_lists  = n(),
    lists    = list(
      tibble(label2 = label2, log2ratio = log2ratio, rank = rank, source = source)
    ),
    .groups = "drop"
  ) %>% 
  filter(regulation == "down") %>% arrange(desc(sum_rank)) %>% 
  head(50) %>% 
  filter(hgnc_symbol %in% c("CYP1B1", "FLNC", "FMO2", "IL1R1", "PRG4")) %>% unnest(lists)
  select(hgnc_symbol, sum_rank) %>% as.data.frame()
  
