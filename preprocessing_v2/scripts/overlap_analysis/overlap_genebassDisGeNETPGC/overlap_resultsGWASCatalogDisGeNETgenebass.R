genebasGrSignatures_overlapChi2
GWASCatalogMoodDisorders_GrSignatures_overlapChi2
disgenetMentalHealth_GrSignatures_overlapChi2

with(
  list(
    A = disgenetMentalHealth_GrSignatures_overlapChi2$processed$original_data$df %>%
      filter(gene_overlap_count > 2,
             p_value < 0.05,
             log2_odds_ratio > 0) %>%
      filter(grepl("global_GR_genes", Var2)) %>%
      pull(overlap_genes) %>%
      strsplit(",") %>%
      unlist() %>%
      str_trim() %>%
      unique(),
    
    B = GWASCatalogMoodDisorders_GrSignatures_overlapChi2$processed$original_data$df %>%
      filter(gene_overlap_count > 2,
             p_value < 0.05,
             log2_odds_ratio > 0) %>%
      filter(grepl("global_GR_genes", Var2)) %>%
      pull(overlap_genes) %>%
      strsplit(",") %>%
      unlist() %>%
      str_trim() %>%
      unique(),
    
    C = genebasGrSignatures_overlapChi2$processed$original_data$df %>%
      filter(gene_overlap_count > 2,
             p_value < 0.05,
             log2_odds_ratio > 0) %>%
      filter(grepl("global_GR_genes", Var2)) %>%
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


GWASCatalog_GrSignatures_overlapChi2$processed$original_data$df %>%
  filter(gene_overlap_count > 2,
         p_value < 0.05,
         log2_odds_ratio > 0) %>%
  filter(grepl("global_GR_genes_", Var2)) %>%
  pull(overlap_genes) %>%
  strsplit(",") %>%
  unlist() %>%
  str_trim() %>%
  unique()



with(
  list(
    A = disgenetMentalHealth_GrSignatures_overlapChi2$processed$original_data$df %>%
      filter(gene_overlap_count > 2,
             p_value < 0.05,
             log2_odds_ratio > 0) %>%
      # filter(grepl("global_GR_genes", Var2)) %>%
      # filter(grepl("BloodCell", Var2)) %>%
      pull(overlap_genes) %>%
      strsplit(",") %>%
      unlist() %>%
      str_trim() %>%
      unique(),
    
    B = GWASCatalogMoodDisorders_GrSignatures_overlapChi2$processed$original_data$df %>%
      filter(gene_overlap_count > 2,
             p_value < 0.05,
             log2_odds_ratio > 0) %>%
      filter(grepl("BloodCell", Var2)) %>%
      pull(overlap_genes) %>%
      strsplit(",") %>%
      unlist() %>%
      str_trim() %>%
      unique(),
    
    C = genebasGrSignatures_overlapChi2$processed$original_data$df %>%
      filter(gene_overlap_count > 2,
             p_value < 0.05,
             log2_odds_ratio > 0) %>%
      # filter(grepl("BloodCell", Var2)) %>%
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
    ABC     = (Reduce(intersect, list(A, B, C)))
  )
)
