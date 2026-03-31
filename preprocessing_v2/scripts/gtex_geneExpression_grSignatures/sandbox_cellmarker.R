grSignatures_cellMarkers2024 = run_enrichr_multi(gene_list = flat_allGrSignatures_18.11.2025[sig_names], database = "CellMarker_2024")

# flat_allGrSignatures_31.10.2025[sig_names]$minusGlobalUpDown5TissuesDerivedCells_BloodCellsUp

grSignatures_cellMarkers2024$enrichr$overlap_only %>% 
  lapply(., head, 10)


grSignatures_cellMarkers2024$enrichr$overlap_only %>%
  imap_dfr(~ .x %>%
             mutate(
               gr_signature = .y,
               tissue = case_when(
                 grepl("liver",  Term, ignore.case = TRUE) ~ "liver",
                 grepl("lung",   Term, ignore.case = TRUE) ~ "lung",
                 grepl("bone",   Term, ignore.case = TRUE) ~ "bone",
                 grepl("kidney", Term, ignore.case = TRUE) ~ "kidney",
                 grepl("gonad",  Term, ignore.case = TRUE) ~ "gonad",
                 grepl("adipose",Term, ignore.case = TRUE) ~ "adipose",
                 grepl("blood",  Term, ignore.case = TRUE) ~ "blood",
                 grepl("brain",  Term, ignore.case = TRUE) ~ "brain",
                 # grepl("Prefrontal Cortex",  Term, ignore.case = TRUE) ~ "brain",
                 grepl("ovary",  Term, ignore.case = TRUE) ~ "ovary",
                 grepl("skin",   Term, ignore.case = TRUE) ~ "skin",
                 TRUE ~ "Other"
               )
             ) %>%
             filter(FDR < 0.05) %>% 
             # filter(tissue != "Other") %>%
             select(gr_signature, Genes, tissue) %>%
             group_by(gr_signature, tissue) %>%
             summarise(
               n_genes = n_distinct(unlist(strsplit(paste(Genes, collapse = ";"), ";"))),
               .groups = "drop"
             )
  ) %>% 
  filter(tissue %in% c("lung", "brain", "blood", "Other")) %>% 
  mutate(gr_signature = str_replace_all(gr_signature,"minusGlobalUpDown5TissuesDerivedCells_", "")) %>% 
  mutate(gr_signature = str_replace_all(gr_signature,"global_GR_genes_globalDown5TissuesDerivedCells", "globalDown")) %>% 
  mutate(gr_signature = str_replace_all(gr_signature,"global_GR_genes_globalUp5TissuesDerivedCells", "globalUp")) %>% 
  as.data.frame() %>% 
  filter(tissue == "brain")


grSignatures_cellMarkers2024$enrichr$overlap_only$minusGlobalUpDown5TissuesDerivedCells_LungCellsDown


grSignatures_cellMarkers2024$enrichr$raw$minusGlobalUpDown5TissuesDerivedCells_BloodCellsUp
 

grSignatures_cellMarkers2024$enrichr$overlap_only %>% lapply(head)


enrichr_minus5$NeuralCells$Up %>% 
  filter(n_genes > 2) %>% 
