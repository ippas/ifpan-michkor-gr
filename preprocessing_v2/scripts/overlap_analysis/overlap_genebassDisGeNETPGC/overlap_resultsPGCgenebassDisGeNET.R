genebasGrSignatures_overlapChi2
pgcGrSignatures_overlapChi2
disgenetMentalHealth_GrSignatures_overlapChi2

with(
  list(
    A = disgenetMentalHealth_GrSignatures_overlapChi2$processed$original_data$df %>%
      filter(gene_overlap_count > 2,
             p_value < 0.05,
             log2_odds_ratio > 0) %>%
      filter(grepl("global_GR_genes_", Var2)) %>%
      pull(overlap_genes) %>%
      strsplit(",") %>%
      unlist() %>%
      str_trim() %>%
      unique(),
    
    B = pgcGrSignatures_overlapChi2$processed$original_data$df %>%
      filter(gene_overlap_count > 2,
             p_value < 0.05,
             log2_odds_ratio > 0) %>%
      filter(Var1 %in% c(
        "pgc_mdd_symptoms_2023-Comm-MDD7_worthless.txt.tsv",
        "pgc_pts_eur_freeze2_overall.results.tsv",
        "pgc_mdd_symptoms_2023-Comm-MDD9_death.txt.tsv",
        "pgc_mdd_symptoms_2023-Comm-MDD1_depressed.txt.tsv",
        "pgc_PGC_MDD2018_10kSNPs.tsv",
        "pgc_mdd_symptoms_2023-Clin-MDD9_death.txt.tsv",
        "pgc_TS_Oct2018.tsv"
      )) %>% 
      filter(grepl("global_GR_genes_", Var2)) %>%
      pull(overlap_genes) %>%
      strsplit(",") %>%
      unlist() %>%
      str_trim() %>%
      unique(),
    
    C = genebasGrSignatures_overlapChi2$processed$original_data$df %>%
      filter(gene_overlap_count > 2,
             p_value < 0.05,
             log2_odds_ratio > 0) %>%
      filter(grepl("global_GR_genes_", Var2)) %>%
      pull(overlap_genes) %>%
      strsplit(",") %>%
      unlist() %>%
      str_trim() %>%
      unique()
  ),
  c(
    A_only  = length(setdiff(A, union(B, C))),
    B_only  = length(setdiff(B, union(A, C))),
    C_only  = length(setdiff(C, union(A, B))),
    AB_only = length(setdiff(intersect(A, B), C)),
    AC_only = length(setdiff(intersect(A, C), B)),
    BC_only = length(setdiff(intersect(B, C), A)),
    ABC     = length(Reduce(intersect, list(A, B, C)))
  )
)


