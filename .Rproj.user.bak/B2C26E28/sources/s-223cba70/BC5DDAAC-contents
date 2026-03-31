

enrichR::listEnrichrDbs() %>% 
  filter(grepl("DisGeNET", libraryName))

enrichR::listEnrichrDbs() %>% 
  filter(grepl("GWAS", libraryName))

enrichR::listEnrichrDbs() %>% 
  filter(grepl("PhenGenI", libraryName))


selected_dbs <- c("DisGeNET", "GWAS_Catalog_2023", "PhenGenI_Association_2021")



genes_DisGeNET <-enrichr_download_gene_lists(
  gene_list = hgnc_symbols_vector_v110,
  database = "DisGeNET")

genes_GWAS_Catalog_2023 <- enrichr_download_gene_lists(
  gene_list = hgnc_symbols_vector_v110,
  database = "GWAS_Catalog_2023")

genes_PhenGenI_Association_2021 <- enrichr_download_gene_lists(
  gene_list = hgnc_symbols_vector_v110,
  database = "PhenGenI_Association_2021")

rbind(
  {genes_DisGeNET %>% 
  filter(n_all_genes > 100) %>% 
  filter(grepl("depression", term, ignore.case = TRUE))},
  {genes_GWAS_Catalog_2023 %>% 
  filter(n_all_genes > 100) %>% 
  filter(grepl("depression", term, ignore.case = TRUE))}) %>% 
  mutate(label = paste(term, Overlap, database, sep = "_")) %>% 
  mutate(label = str_replace_all(label, " ", "_")) %>% 
  ungroup %>% 
  as.data.frame() %>% 
  select(label, combined_Genes)

genes_PhenGenI_Association_2021 %>% 
  filter(n_all_genes > 50) %>% 
  filter(grepl("depression", Term, ignore.case = TRUE))





