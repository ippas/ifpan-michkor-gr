flat_allGrSignatures_17.10.2025$global_GR_genes_globalUp5TissuesDerivedCells

dbs <- listEnrichrDbs()
dbs %>% 
  filter(grepl("cell", libraryName, ignore.case = T))

tmp_enrichr <- run_enrichr(gene_list = flat_allGrSignatures_17.10.2025$global_GR_genes_globalUp5TissuesDerivedCells,
            database = "ChEA_2022")

tmp_enrichr %>% 
  filter(grepl("NR3C1", Term, ignore.case = T)) %>% 
  .$Genes %>% 
  strsplit(";") %>% 
  unlist %>% unique() %>% length()


tmp_enrichr <- run_enrichr(gene_list = flat_allGrSignatures_17.10.2025$global_GR_genes_globalUp4TissuesDerivedCells,
                           database = "ChEA_2022")

tmp_enrichr %>% 
  filter(grepl("NR3C1", Term, ignore.case = T)) %>% 
  .$Genes %>% 
  strsplit(";") %>% 
  unlist %>% unique() %>% length()

tmp_enrichr <- run_enrichr(gene_list = flat_allGrSignatures_17.10.2025$minusGlobalUp5TissuesDerivedCells_LungCellsUp,
                           database = "CellMarker_2024")

tmp_enrichr %>% 
  filter(grepl("Lung|Blood|Brain", Term, ignore.case = T)) %>% 
  filter(Adjusted.P.value < 0.05) %>% 
  filter(n_genes > 2)


tmp_enrichr <- run_enrichr(gene_list = flat_allGrSignatures_17.10.2025$minusGlobalUp5TissuesDerivedCells_NeuralCellsUp,
                           database = "CellMarker_2024")

tmp_enrichr %>% 
  filter(grepl("Lung|Blood|Brain", Term, ignore.case = T)) %>% 
  filter(Adjusted.P.value < 0.05) %>% 
  filter(n_genes > 2)


tmp_enrichr <- run_enrichr(gene_list = c(flat_allGrSignatures_17.10.2025$minusGlobalUp5TissuesDerivedCells_LungCellsUp, 
                                         flat_allGrSignatures_17.10.2025$minusGlobalDown5TissuesDerivedCells_LungCellsDown),
                           database = "CellMarker_2024")

tmp_enrichr %>% 
  filter(Adjusted.P.value < 0.05) %>% 
  filter(n_genes > 2)  %>% 
  arrange(desc(Combined.Score)) %>% 
  filter(grepl("Blood", Term)) %>%  
  .$Genes %>% strsplit(";") %>% unlist %>% unique



flat_allGrSignatures_17.10.2025$global_GR_genes_globalUp5TissuesDerivedCells


tmp_enrichr %>% 
  filter(n_genes > 2) %>% 
  filter(Adjusted.P.value < 0.05)



# ##############################################################################
flat_allGrSignatures_17.10.2025$minusGlobalUp4TissuesDerivedCells_LungCellsUp %>% length()


c(flat_allGrSignatures_17.10.2025$minusClusterD_NeuralCellsDown,
  flat_allGrSignatures_17.10.2025$minusClusterD_LungCellsDown,
  flat_allGrSignatures_17.10.2025$minusClusterD_BloodCellsDown ) %>%  
  table %>% as.data.frame() %>% filter(Freq > 1)
