# prepare gene lists form MGI, using protein coding genes from biomart v110

enrichR::listEnrichrDbs() %>% 
  filter(grepl("MGI", libraryName))


genes_MGI_2024 <- enrichr_download_gene_lists(
  gene_list = hgnc_symbols_vector_v110,
  database = "MGI_Mammalian_Phenotype_Level_4_2024"
)

genes_MGI_2024 %>% 
  mutate(genes = str_replace_all(combined_genes, "; ", "|")) %>% 
  select(-c(combined_genes, database)) -> genes_MGI_2024


genes_MGI_2024 %>% 
  filter(grepl("lipid|cholesterol|LDL|HDL", term, ignore.case = T)) %>% 
  .$genes %>% 
  str_split("\\|") %>% 
  unlist %>% 
  unique() -> MGI_lipid_gene_signature


genes_MGI_2024 %>% 
  filter(grepl("glucose", term, ignore.case = T)) %>% 
  .$genes %>% 
  str_split("\\|") %>% 
  unlist %>% 
  unique() -> MGI_glucose_gene_signature

overlap_coefficient(set1 = MGI_lipid_gene_signature, set2 = MGI_glucose_gene_signature)
