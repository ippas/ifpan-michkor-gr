
# wczytanie
read_tsv("data/databases/pgc-genes-gwas-combine-filtered.tsv", col_names = T) %>% 
  mutate(label = str_replace(file, ".meta.gz", "")) %>% 
  mutate(label = paste0(label, "_", pvalue_threshold)) %>% 
  select(c(label, unique_genes))



# utworzenie list 
read_tsv("data/databases/pgc-genes-gwas-combine-filtered.tsv", col_names = TRUE) %>%
  # Tworzenie etykiety
  mutate(label = str_replace(file, "\\.meta\\.gz$", "")) %>% 
  mutate(label = paste0(label, "_", pvalue_threshold)) %>%
  # Przetworzenie unique_genes na listę genów
  mutate(unique_genes = str_split(unique_genes, ",")) %>%  # Zakładając, że geny są oddzielone przecinkami
  # Grupowanie po label i agregacja genów
  group_by(label) %>%
  summarise(unique_genes = list(unique(unlist(unique_genes)))) %>%
  ungroup() %>%
  # Konwersja do listy z nazwami
  deframe() -> pgc_gene_lists




# overlap