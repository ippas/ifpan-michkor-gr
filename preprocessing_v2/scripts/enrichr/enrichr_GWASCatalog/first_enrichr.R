# ##############################################################################
# ---- uses data ----
# ##############################################################################

AllBrain2BrainSignatures2Global


# ##############################################################################
# ---- run analysis ----
# ##############################################################################
AllBrain2BrainSignatures2Global_GWAS_Enrichr <- 
  AllBrain2BrainSignatures2Global %>% 
  lapply(function(x) {
    run_enrichr(x, database = "GWAS_Catalog_2023")
  })


AllBrain2BrainSignatures2Global_GWAS_Enrichr %>% 
  lapply(., function(x){
    x %>% 
      separate(Overlap, into = c("n_genes", "n_total"), sep = "/", convert = TRUE) %>%
      mutate(n_genes = as.integer(n_genes)) %>% 
      filter(P.value < 0.05) %>% 
      filter(n_genes > 2)
  })
