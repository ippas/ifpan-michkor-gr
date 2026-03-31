summarize_pgc_overlap <- function(df,
                                  exclude_var1 = c("brain_up", "brain_down", 
                                                   "metasignature_up", "metasignature_down"),
                                  p_value_threshold = 0.05,
                                  min_gene_overlap = 3) {
  library(dplyr)
  
  # jeśli kilka progów -> iteracja
  if (length(p_value_threshold) > 1) {
    results <- lapply(p_value_threshold, function(thr) {
      summarize_pgc_overlap(df = df,
                            exclude_var1 = exclude_var1,
                            p_value_threshold = thr,
                            min_gene_overlap = min_gene_overlap)
    })
    names(results) <- as.character(p_value_threshold)
    
    # zbiorcza tabela
    summary_df <- do.call(rbind, lapply(results, function(x) x$summary_row))
    rownames(summary_df) <- NULL
    
    return(list(
      results_per_threshold = results,
      summary_table = summary_df
    ))
  }
  
  # ------------------------
  # filtrowanie
  # ------------------------
  df_filtered <- df %>%
    filter(!(Var1 %in% exclude_var1)) %>%
    filter(gene_overlap_count >= min_gene_overlap) %>%
    filter(p_value < p_value_threshold) %>%
    select(-fdr)
  
  # wszystkie potencjalne Var1 / Var2
  all_Var1 <- unique(df$Var1[!(df$Var1 %in% exclude_var1)])
  all_Var2 <- unique(df$Var2)
  
  # które pozostały po filtrach
  kept_Var1 <- unique(df_filtered$Var1)
  kept_Var2 <- unique(df_filtered$Var2)
  
  # które całkowicie wypadły
  excluded_Var1 <- setdiff(all_Var1, kept_Var1)
  excluded_Var2 <- setdiff(all_Var2, kept_Var2)
  
  # ------------------------
  # liczby podstawowe
  # ------------------------
  n_association <- nrow(df_filtered)
  n_genes <- df_filtered$overlap_genes %>% strsplit(",") %>% unlist() %>% unique() %>% length()
  n_papers <- df_filtered %>%
    select(Var1) %>% unique() %>%
    mutate(source = sapply(strsplit(Var1, "_"), `[`, 1)) %>%
    select(source) %>% unique() %>% nrow()
  n_phenotypes <- length(unique(df_filtered$Var2))
  n_geneLists <- length(unique(df_filtered$Var1))
  
  # fenotypy
  pheno_counts <- df_filtered %>%
    count(Var2, name = "Freq") %>%
    arrange(desc(Freq))
  n_single_assoc_pheno <- sum(pheno_counts$Freq == 1)
  top_pheno <- pheno_counts %>% filter(Freq == max(Freq)) %>% pull(Var2)
  top_pheno_genes <- df_filtered %>%
    filter(Var2 %in% top_pheno) %>%
    .$overlap_genes %>% strsplit(",") %>% unlist() %>% unique()
  top_pheno_nGenes <- length(top_pheno_genes)
  
  # geneList
  gl_counts <- df_filtered %>%
    count(Var1, name = "Freq") %>%
    arrange(desc(Freq))
  top_geneList <- gl_counts %>% filter(Freq == max(Freq)) %>% pull(Var1)
  top_geneList_genes <- df_filtered %>%
    filter(Var1 %in% top_geneList) %>%
    .$overlap_genes %>% strsplit(",") %>% unlist() %>% unique()
  top_geneList_nGenes <- length(top_geneList_genes)
  
  # per source
  df_filtered_source <- df_filtered %>%
    mutate(source = sapply(strsplit(Var1, "_"), `[`, 1))
  pheno_per_source <- df_filtered_source %>%
    select(Var2, source) %>% unique()
  pheno_source_counts <- pheno_per_source %>% count(Var2, name = "Freq")
  top_pheno_source <- pheno_source_counts %>% filter(Freq == max(Freq)) %>% pull(Var2)
  n_sources_for_top_pheno <- ifelse(nrow(pheno_source_counts) > 0, max(pheno_source_counts$Freq), 0)
  
  # ------------------------
  # podsumowanie jako jeden wiersz
  # ------------------------
  summary_row <- data.frame(
    pval_threshold = p_value_threshold,
    n_association = n_association,
    n_genes = n_genes,
    n_papers = n_papers,
    n_phenotypes = n_phenotypes,
    n_geneLists = n_geneLists,
    n_single_assoc_pheno = n_single_assoc_pheno,
    top_pheno = paste(top_pheno, collapse = ";"),
    top_pheno_nGenes = top_pheno_nGenes,
    top_geneList = paste(top_geneList, collapse = ";"),
    top_geneList_nGenes = top_geneList_nGenes,
    top_pheno_per_source = paste(top_pheno_source, collapse = ";"),
    n_sources_for_top_pheno = n_sources_for_top_pheno,
    stringsAsFactors = FALSE
  )
  
  # ------------------------
  # zwrot
  # ------------------------
  return(list(
    summary_row = summary_row,
    filtered_df = df_filtered,
    pheno_counts = pheno_counts,
    gl_counts = gl_counts,
    pheno_source_counts = pheno_source_counts,
    excluded_Var1 = excluded_Var1,
    excluded_Var2 = excluded_Var2
  ))
}


list(
  res1e4 = summarize_pgc_overlap(
    pgcGeneCenter50kb_AllBrain2Brain2Global_rsIDp1e4$overlap$original_data$df,
    p_value_threshold = c(0.05, 0.01, 0.001, 1e-4, 1e-5),
    min_gene_overlap = 3
  ),
  res1e5 = summarize_pgc_overlap(
    pgcGeneCenter50kb_AllBrain2Brain2Global_rsIDp1e5$overlap$original_data$df,
    p_value_threshold = c(0.05, 0.01, 0.001, 1e-4, 1e-5),
    min_gene_overlap = 3
  ),
  res1e6 = summarize_pgc_overlap(
    pgcGeneCenter50kb_AllBrain2Brain2Global_rsIDp1e6$overlap$original_data$df,
    p_value_threshold = c(0.05, 0.01, 0.001, 1e-4, 1e-5),
    min_gene_overlap = 3
  ),
  res1e7 = summarize_pgc_overlap(
    pgcGeneCenter50kb_AllBrain2Brain2Global_rsIDp1e7$overlap$original_data$df,
    p_value_threshold = c(0.05, 0.01, 0.001, 1e-4, 1e-5),
    min_gene_overlap = 3
  ),
  res1e8 = summarize_pgc_overlap(
    pgcGeneCenter50kb_AllBrain2Brain2Global_rsIDp1e8$overlap$original_data$df,
    p_value_threshold = c(0.05, 0.01, 0.001, 1e-4, 1e-5),
    min_gene_overlap = 3
  )
) -> tmp

imap_dfr(tmp, ~ {
  .x$summary_table %>%
    mutate(source = .y)   # dodaje nazwę listy (np. "res1e4")
}) %>% 
  select(c(source, n_genes, pval_threshold, n_association, n_papers, n_phenotypes, n_geneLists, n_single_assoc_pheno)) %>% 
  ggplot(aes(x = pval_threshold, y = n_genes, 
             color = source, group = source)) +
  geom_point(size = 3) +
  geom_line() +
  scale_x_log10(breaks = c(1e-2, 1e-3, 1e-4, 1e-5)) +   # log-skala, czytelniejsza dla pval
  theme_minimal(base_size = 14) +
  labs(
    x = "P-value threshold",
    y = "Number of genes",
    color = "Source",
    title = "Liczba genów w zależności od progu istotności"
  )




tmp$res1e5$results_per_threshold$`0.001`$filtered_df
