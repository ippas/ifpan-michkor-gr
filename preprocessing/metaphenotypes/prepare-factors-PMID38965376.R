
file <- "data/metaphenotypes/PMID38965376_FA_STable5_41562_2024_1909_MOESM5_ESM.xlsx"

latentFactors_pmid38965376 <- setNames(
  lapply(excel_sheets(file)[-1], \(s) read_excel(file, sheet = s)),
  excel_sheets(file)[-1]
)

latentFactors_pmid38965376_df <- bind_rows(latentFactors_pmid38965376, .id = "latent_factor") %>% 
  set_colnames(c("latent_factor", "item", "item_name", "est", "se", "z", "p"))

latentFactors_pmid38965376_df %>% 
  group_by(latent_factor) %>% 
  nest %>% 
  mutate(n_phenotypes = map(data, ~nrow(.x))) %>% 
  unnest(n_phenotypes)


latentFactors_pmid38965376_df %>% head

latentFactors_pmid38965376_df$item_name %>% 
  unique %>% length()

################################################################################
file <- "data/metaphenotypes/simpleEnrichr_grSignatures_GWASCatalog.xlsx"

enrichr_grSignatures_GWASCatalog <- setNames(
  lapply(excel_sheets(file), \(s) read_excel(file, sheet = s)),
  excel_sheets(file)
)

enrichr_grSignatures_GWASCatalog[1:4] %>% 
  bind_rows(.id = "grSignature") -> enrichr_grSignatures_GWASCatalog_df

enrichr_grSignatures_GWASCatalog_df$Phenotype %>% unique() %>% length()



intersect(
  tolower(unique(latentFactors_pmid38965376_df$item_name)),
  tolower(unique(enrichr_grSignatures_GWASCatalog_df$Phenotype))
)
