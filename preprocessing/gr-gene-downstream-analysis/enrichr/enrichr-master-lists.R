install.packages("enrichR")
require(enrichR)

dbs <- listEnrichrDbs()

enrichr_pathway_names <- c("BioPlanet_2019", "KEGG_2021_Human", "WikiPathway_2023_Human", "Elsevier_Pathway_Collection")
enrichr_diseases_drugs_names <- c("DSigDB", "DrugMatrix", "GeDiPNet_2023", "DisGeNET")
enrichr_ontologies <- c("GO_Biological_Process_2023")

gr_master_up

dbs %>% 
  filter(libraryName %in% c(enrichr_pathway_names, enrichr_diseases_drugs_names, enrichr_ontologies)) -> dbs_filtered

dbs_filtered %>% dim

lapply(c(1:9), function(i){
  enrichr(gr_master_up$hgnc_symbol, dbs_filtered[i,])
}) -> enriched_gr_master_up

names(enriched_gr_master_up) <- c(enrichr_pathway_names, enrichr_diseases_drugs_names, enrichr_ontologies)

lapply(c(1:9), function(i){
  enriched_gr_master_up[i][[1]][3]
}) %>% unlist(recursive = F) -> enriched_gr_master_up

enriched_gr_master_up %>% 
  bind_rows(., .id = "enrichr_database") %>% 
  as.tibble() %>% 
  mutate(n_genes = map(Genes, ~length(convert_genes_to_vector(.x, split = ";")))) %>% 
  unnest(n_genes) %>% as.data.frame() %>% 
  filter(n_genes >= 3) %>%
  filter(Adjusted.P.value < 0.01) %>% 
  filter(enrichr_database %in% c("DSigDB"))

enriched_gr_master_up %>% 
  bind_rows(., .id = "enrichr_database") %>% 
  as.tibble() %>% 
  mutate(n_genes = map(Genes, ~length(convert_genes_to_vector(.x, split = ";")))) %>% 
  unnest(n_genes) %>% as.data.frame() %>% 
  filter(n_genes >= 3) %>%
  filter(Adjusted.P.value < 0.01) %>%  
  # filter(enrichr_database %in% enrichr_diseases_drugs_names) %>% 
  write.table(., file = "results/google-drive/enrichr-master-up-top50-ranked-score-fdr0.01-overlap3.tsv", sep = "\t", col.names = T, row.names = F, quote = F)
  
  
