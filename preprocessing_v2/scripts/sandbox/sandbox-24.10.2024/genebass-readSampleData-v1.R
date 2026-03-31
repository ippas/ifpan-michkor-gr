

genebassData_SKATO_p0.0001 <- read.delim("data/genebass_v2/merged_and_filtered_genebassData/genebass_allCategories_SKATO_mergedFilteredP0001.tsv.bgz")


genebassData_SKATO_p0.0001 %>% 
  select(-heritability, -pvalue_threshold, -pvalue_test) %>% 
  filter(annotation == "pLoF") %>% 
  arrange(phenocode) %>% 
  group_by(phenocode) %>% nest %>% 
  mutate(n_genes = map(data, ~ .x$gene_symbol %>% unique %>% length)) %>% 
  unnest(n_genes) %>% 
  filter(n_genes > 10) %>% 
  unnest(data) %>% 
  as.data.frame() %>% 
  # mutate(description  = str_replace_all(" ", "_", description )) %>% 
  filter(gene_symbol %in% flat_allGrSignatures_17.10.2025$NeuralCellsUp) %>% 
  filter(description == "Polyunsaturated fat")
  
  
genebassData_SKATO_p0.0001 %>% 
  filter(grepl("mental", description)) %>% 
  filter(description == "Polyunsaturated fat")
