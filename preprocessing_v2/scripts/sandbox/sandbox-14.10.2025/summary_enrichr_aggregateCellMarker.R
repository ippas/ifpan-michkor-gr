# ##############################################################################
# ---- uses data ----
# ##############################################################################

gene_summaryList_blood

gene_summaryList_brain

gene_summaryList_lung$n10_short_time$summary_up$enrichr_summary_source$enrichr_3pub


# ##############################################################################
gene_summaryList_blood$n10_short_time$summary_up$enrichr_summary_source$genes_3pub %>% length()

gene_summaryList_blood$n10_short_time$summary_down$enrichr_summary_source$enrichr_3pub$CellMarker_2024 %>%
  filter(Adjusted.P.value < 0.05) %>%
  filter(n_genes >= 3) %>%
  filter(grepl("lung|blood|brain", Term, ignore.case = TRUE)) %>%
  mutate(
    aggregate_term = case_when(
      grepl("lung", Term, ignore.case = TRUE)  ~ "lung",
      grepl("blood", Term, ignore.case = TRUE) ~ "blood",
      grepl("brain", Term, ignore.case = TRUE) ~ "brain",
      TRUE ~ "other"
    )
  ) %>% 
  group_by(aggregate_term) %>% 
  nest %>% 
  mutate(aggregate_genes = map(data, ~ .x %>% .$Genes %>% unlist %>% strsplit(";") %>% unlist %>% unique)) %>% 
  mutate(n_aggreate_genes = map(aggregate_genes, ~ .x %>% length)) %>% 
  unnest(n_aggreate_genes) %>% 
  mutate(n_traits = map(data, ~ .x %>% nrow)) %>% 
  unnest(n_traits)


gene_summaryList_blood$n10_short_time$summary_up$enrichr_summary_source$enrichr_3pub$ChEA_2022 %>%
  filter(Adjusted.P.value < 0.05) %>%
  filter(n_genes >= 3) %>%
  filter(grepl("NR3C1|NR3C2", Term, ignore.case = TRUE)) %>% 
  mutate(
    aggregate_term = case_when(
      grepl("NR3C1", Term, ignore.case = TRUE)  ~ "NR3C1",
      grepl("NR3C2", Term, ignore.case = TRUE) ~ "NR3C2",
      TRUE ~ "other"
    )
  ) %>% 
  group_by(aggregate_term) %>% 
  nest %>% 
  mutate(aggregate_genes = map(data, ~ .x %>% .$Genes %>% unlist %>% strsplit(";") %>% unlist %>% unique)) %>% 
  mutate(n_aggreate_genes = map(aggregate_genes, ~ .x %>% length)) %>% 
  unnest(n_aggreate_genes) %>% 
  mutate(n_traits = map(data, ~ .x %>% nrow)) %>% 
  unnest(n_traits)


gene_summaryList_blood$n10_short_time$summary_down$enrichr_summary_source$enrichr_3pub$LINCS_L1000_Chem_Pert_Consensus_Sigs %>%
  filter(Adjusted.P.value < 0.05) %>%
  filter(n_genes >= 3) %>%  
  filter(grepl("Betamethasone|Dexamethasone|Fluocortolone|Methylprednisolone|Paramethasone|Prednisolone|Prednisone|Triamcinolone|Hydrocortisone|Cortisone|Prednylidene|Rimexolone|Deflazacort|Cloprednol|Meprednisone|Cortivazol|Vamorolone", Term, ignore.case = TRUE)) %>%
  mutate(
    aggregate_term = case_when(
      grepl("up", Term, ignore.case = TRUE)  ~ "GC_UP_NIPHID_H02AB",
      grepl("down", Term, ignore.case = TRUE) ~ "GC_DOWN_NIPHID_H02AB",
      TRUE ~ "other"
    )
  ) %>% 
  group_by(aggregate_term) %>% 
  group_by(aggregate_term) %>% 
  nest %>% 
  mutate(aggregate_genes = map(data, ~ .x %>% .$Genes %>% unlist %>% strsplit(";") %>% unlist %>% unique)) %>% 
  mutate(n_aggreate_genes = map(aggregate_genes, ~ .x %>% length)) %>% 
  unnest(n_aggreate_genes) %>% 
  mutate(n_traits = map(data, ~ .x %>% nrow)) %>% 
  unnest(n_traits)
