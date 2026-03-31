# ##############################################################################
# ---- functions ----
# ##############################################################################
compute_disgenet_pgc_genebass_venn <- function(
    pattern_var2 = NULL    # <- nowy argument do filtra po Var2
) {
  
  # ------------------------------
  # 1) DisGeNET
  # ------------------------------
  df_dis <- disgenetMentalHealth_GrSignatures_overlapChi2$processed$original_data$df
  
  if (!is.null(pattern_var2)) {
    df_dis <- df_dis %>% filter(grepl(pattern_var2, Var2))
  }
  
  genes_disgenet <- df_dis %>% 
    filter(gene_overlap_count > 2, p_value < 0.05) %>% 
    pull(overlap_genes) %>% 
    strsplit(",") %>% unlist() %>% unique()
  
  
  # ------------------------------
  # 2) PGC
  # ------------------------------
  df_pgc <- pgcGrSignatures_overlapChi2$processed$original_data$df
  
  if (!is.null(pattern_var2)) {
    df_pgc <- df_pgc %>% filter(grepl(pattern_var2, Var2))
  }
  
  genes_pgc <- df_pgc %>%
    filter(gene_overlap_count > 2, p_value < 0.05) %>% 
    filter(Var1 %in% c(
      "pgc_mdd_symptoms_2023-Comm-MDD7_worthless.txt.tsv",
      "pgc_pts_eur_freeze2_overall.results.tsv",
      "pgc_mdd_symptoms_2023-Comm-MDD9_death.txt.tsv",
      "pgc_mdd_symptoms_2023-Comm-MDD1_depressed.txt.tsv",
      "pgc_PGC_MDD2018_10kSNPs.tsv",
      "pgc_mdd_symptoms_2023-Clin-MDD9_death.txt.tsv",
      "pgc_TS_Oct2018.tsv"
    )) %>%
    pull(overlap_genes) %>%
    strsplit(",") %>% unlist() %>% unique()
  
  
  # ------------------------------
  # 3) GeneBass
  # ------------------------------
  df_gb <- genebasGrSignatures_overlapChi2$processed$original_data$df
  
  if (!is.null(pattern_var2)) {
    df_gb <- df_gb %>% filter(grepl(pattern_var2, Var2))
  }
  
  genes_genebass <- df_gb %>% 
    filter(gene_overlap_count > 2, p_value < 0.05) %>% 
    pull(overlap_genes) %>%
    strsplit(",") %>% unlist() %>% unique() %>% sort()
  
  
  # ------------------------------
  # 4) Obszary Venn
  # ------------------------------
  all3 <- Reduce(intersect, list(genes_disgenet, genes_pgc, genes_genebass))
  
  Dis_PGC_only  <- setdiff(intersect(genes_disgenet, genes_pgc),     all3)
  Dis_GB_only   <- setdiff(intersect(genes_disgenet, genes_genebass),all3)
  PGC_GB_only   <- setdiff(intersect(genes_pgc,      genes_genebass),all3)
  
  Dis_only      <- setdiff(genes_disgenet, union(genes_pgc, genes_genebass))
  PGC_only      <- setdiff(genes_pgc,      union(genes_disgenet, genes_genebass))
  GB_only       <- setdiff(genes_genebass, union(genes_disgenet, genes_pgc))
  
  
  # ------------------------------
  # 5) Listy genów
  # ------------------------------
  genes_list <- list(
    DisGeNET_only          = Dis_only,
    PGC_only               = PGC_only,
    GeneBass_only          = GB_only,
    DisGeNET_PGC_only      = Dis_PGC_only,
    DisGeNET_GeneBass_only = Dis_GB_only,
    PGC_GeneBass_only      = PGC_GB_only,
    All3                   = all3
  )
  
  # ------------------------------
  # 6) Liczniki
  # ------------------------------
  counts_list <- lapply(genes_list, length)
  
  # ------------------------------
  # 7) Zwracamy
  # ------------------------------
  list(
    genes  = genes_list,
    counts = counts_list
  )
}

compute_disgenet_pgc_genebass_venn(pattern_var2 = NULL)

compute_disgenet_pgc_genebass_venn(pattern_var2 = "minusGlobalUpDown5TissuesDerivedCells_NeuralCellsUp|minusGlobalUpDown5TissuesDerivedCells_NeuralCellsDown")


disgenetMentalHealth_GrSignatures_overlapChi2$processed$original_data$df$Var2 %>% unique

# # ##############################################################################
# 
# # ---- 1. Wyciąganie listy genów z DisGeNET ----
# genes_disgenet <- disgenetMentalHealth_GrSignatures_overlapChi2$processed$original_data$df %>% 
#   filter(gene_overlap_count > 2, p_value < 0.05) %>% 
#   pull(overlap_genes) %>% 
#   strsplit(",") %>% 
#   unlist() %>% 
#   unique()
# 
# 
# genes_pgc <- pgcGrSignatures_overlapChi2$processed$original_data$df %>% 
#   filter(gene_overlap_count > 2, p_value < 0.05) %>% 
#   filter(Var1 %in% c("pgc_mdd_symptoms_2023-Comm-MDD7_worthless.txt.tsv",
#                      "pgc_pts_eur_freeze2_overall.results.tsv",
#                      "pgc_mdd_symptoms_2023-Comm-MDD9_death.txt.tsv",
#                      "pgc_mdd_symptoms_2023-Comm-MDD1_depressed.txt.tsv",
#                      "pgc_PGC_MDD2018_10kSNPs.tsv",
#                      "pgc_mdd_symptoms_2023-Clin-MDD9_death.txt.tsv",
#                      "pgc_TS_Oct2018.tsv"
#                      # "pgc_PGC3_SCZ_wave3.asian.autosome.public.v3.vcf.tsv.tsv",
# 
#                      )) %>% 
#   .$overlap_genes %>% 
#   strsplit(",") %>% unlist %>% unique
# 
# 
# genebasGrSignatures_overlapChi2$processed$original_data$df %>% 
#   filter(gene_overlap_count > 2) %>% 
#   filter(p_value < 0.05) %>% 
#   .$overlap_genes %>% 
#   strsplit(",") %>% unlist %>% 
#   unique %>% sort -> genes_genebass
# 
# 
# 
# 
# 
# # ---- 2. Wyciąganie listy genów z PGC ----
# # genes_pgc <- pgcGrSignatures_overlapChi2$processed$original_data$df %>% 
# #   filter(gene_overlap_count > 2, p_value < 0.05) %>% 
# #   pull(overlap_genes) %>% 
# #   strsplit(",") %>% 
# #   unlist() %>% 
# #   unique()


disgenetMentalHealth_GrSignatures_overlapChi2$processed$original_data$df %>% 
  filter(gene_overlap_count > 2) %>% 
  filter(p_value < 0.05) %>% 
  filter(grepl("Schizophrenia", Var1)) %>% 
  filter(gene_overlap_count == 16)
