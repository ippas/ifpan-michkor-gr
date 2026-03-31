res <- run_enrichr_multi(
  gene_lists = flat_allGrSignatures_31.10.2025[sig_names],
  database = "GWAS_Catalog_2025",
  min_overlap_genes = 3,
  fdr_threshold = 0.05,
  names_mapper = names_mapper,
  xlsx_file = "results_v2/overlap/enrichr/enrichr_GWASCatalog2025_GRsignatures_11.01.2026.xlsx"
)

res$enrichr$raw %>% 
  lapply(., function(x){
    x %>% filter(n_genes > 2) %>% 
      # filter(FDR < 0.01)
      filter(n_genes > 10)
      # head(10)
  }) %>% 
  lapply(., function(x){
    x$Term
  }) %>% unname() %>% unlist %>% unique -> term_GWASCatalog2025


res$enrichr$raw %>% 
  lapply(., function(x){
    x %>% filter(Term %in% term_GWASCatalog2025)  
  }) 
