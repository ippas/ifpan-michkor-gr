# ##############################################################################
# ---- ueses data ----
# ##############################################################################
grSignature_genebass_associations_FDRMonteCarlo$original_genebass_association$brain_up

genebass_mentalHealth_skat_all$annotation %>% unique

# ##############################################################################
# ---- summary results ----
# ##############################################################################
signatures <- c("brain_up", "brain_down", "metasignature_up", "metasignature_down")

results_df <- lapply(signatures, function(sig_name) {
  df <- grSignature_genebass_associations_FDRMonteCarlo$original_genebass_association[[sig_name]] %>%
    select(!c(p40_pvalue, p60_pvalue, p80_pvalue, p90_pvalue, q1_pvalue, q3_pvalue, p30_pvalue)) %>%
    filter(pvalue < 0.05) %>%
    mutate(
      fdr_0.2 = pvalue < p20_pvalue,
      fdr_0.1 = pvalue < p10_pvalue
    )
  
  df_fdr01 <- df %>% filter(fdr_0.1)
  df_p1e5  <- df %>% filter(pvalue < 1e-5 & fdr_0.1)
  
  data.frame(
    signature = sig_name,
    
    # FDR < 0.1
    n_variants = nrow(df_fdr01),
    n_genes = length(unique(df_fdr01$gene_symbol)),
    n_phenotypes = length(unique(df_fdr01$phenocode)),
    n_categories = length(unique(df_fdr01$category)),
    
    # p < 1e-5 AND FDR < 0.1
    n_genes_p1e5 = length(unique(df_p1e5$gene_symbol)),
    n_variants_p1e5 = nrow(df_p1e5),
    n_phenotypes_p1e5 = length(unique(df_p1e5$phenocode)),
    genes_p1e5 = paste(sort(unique(df_p1e5$gene_symbol)), collapse = "; ")
  )
}) %>% bind_rows()
results_df



results_df <- lapply(signatures, function(sig_name) {
  df <- grSignature_genebass_associations_FDRMonteCarlo$original_genebass_association[[sig_name]] %>%
    select(!c(p40_pvalue, p60_pvalue, p80_pvalue, p90_pvalue, q1_pvalue, q3_pvalue, p30_pvalue)) %>%
    filter(pvalue < 0.05) %>%
    mutate(
      fdr_0.2 = pvalue < p20_pvalue,
      fdr_0.1 = pvalue < p10_pvalue
    ) %>%
    filter(fdr_0.1)
  
  # rozdziel według kategorii
  df %>%
    group_by(category) %>%
    group_split() %>%
    lapply(function(df_cat) {
      df_p1e5 <- df_cat %>% filter(pvalue < 1e-5)
      
      data.frame(
        signature = sig_name,
        category = unique(df_cat$category),
        
        n_variants = nrow(df_cat),
        n_genes = length(unique(df_cat$gene_symbol)),
        n_phenotypes = length(unique(df_cat$phenocode)),
        
        n_genes_p1e5 = length(unique(df_p1e5$gene_symbol)),
        n_variants_p1e5 = nrow(df_p1e5),
        n_phenotypes_p1e5 = length(unique(df_p1e5$phenocode)),
        genes_p1e5 = paste(sort(unique(df_p1e5$gene_symbol)), collapse = "; ")
      )
    }) %>% bind_rows()
}) %>% bind_rows()

results_df %>% 
  filter(n_genes_p1e5 > 0)


data.frame(
  n_variants = lapply(c("brain_up", "brain_down", "metasignature_up", "metasignature_down"), function(sig)
    grSignature_genebass_associations_FDRMonteCarlo$original_genebass_association[[sig]] %>%
      filter(pvalue < 0.05, pvalue < p10_pvalue)
  ) %>% bind_rows() %>% nrow(),
  
  n_genes = lapply(c("brain_up", "brain_down", "metasignature_up", "metasignature_down"), function(sig)
    grSignature_genebass_associations_FDRMonteCarlo$original_genebass_association[[sig]] %>%
      filter(pvalue < 0.05, pvalue < p10_pvalue) %>%
      pull(gene_symbol)
  ) %>% unlist() %>% unique() %>% length(),
  
  n_phenotypes = lapply(c("brain_up", "brain_down", "metasignature_up", "metasignature_down"), function(sig)
    grSignature_genebass_associations_FDRMonteCarlo$original_genebass_association[[sig]] %>%
      filter(pvalue < 0.05, pvalue < p10_pvalue) %>%
      pull(phenocode)
  ) %>% unlist() %>% unique() %>% length()
)


# ##############################################################################
# 1) build one big table with signature & random_fdr (pvalue < 0.05)
all_assoc <- map_df(signatures, function(sig_name) {
  grSignature_genebass_associations_FDRMonteCarlo$original_genebass_association[[sig_name]] %>%
    filter(pvalue < 0.05) %>%
    select(-c(q1_pvalue, q3_pvalue, mean_pvalue, p30_pvalue, p40_pvalue, p60_pvalue, p80_pvalue, p90_pvalue)) %>%
    mutate(
      signature  = sig_name,
      random_fdr = case_when(
        pvalue < p5_pvalue  ~ "0.05",
        pvalue < p10_pvalue ~ "0.1",
        pvalue < p20_pvalue ~ "0.2",
        TRUE                ~ "ns"
      )
    )
})

# 1) Zbuduj all_assoc (wszystkie sygnatury, p < 0.05) z random_fdr
all_assoc <- map_df(signatures, function(sig_name) {
  grSignature_genebass_associations_FDRMonteCarlo$original_genebass_association[[sig_name]] %>%
    filter(pvalue < 0.05) %>%
    select(-c(q1_pvalue, q3_pvalue, mean_pvalue, p30_pvalue, p40_pvalue, p60_pvalue, p80_pvalue, p90_pvalue)) %>%
    mutate(
      signature  = sig_name,
      random_fdr = case_when(
        pvalue < p5_pvalue  ~ "0.05",
        pvalue < p10_pvalue ~ "0.1",
        pvalue < p20_pvalue ~ "0.2",
        TRUE                ~ "ns"
      )
    )
})

# 2) Podsumowanie dla p < 0.05 (broad)
summary_broad <- all_assoc %>%
  group_by(signature, random_fdr) %>%
  summarise(
    n_variants_broad    = length(signature),
    n_genes_broad       = n_distinct(gene_symbol),
    n_phenotypes_broad  = n_distinct(phenocode),
    n_categories_broad  = n_distinct(category),
    .groups             = "drop"
  )

# 3) Podsumowanie dla p < 1e-5 (strict), łącznie z listą genów
summary_strict <- all_assoc %>%
  filter(pvalue < 1e-5) %>%
  group_by(signature, random_fdr) %>%
  summarise(
    n_variants_strict   = length(signature),
    n_genes_strict      = n_distinct(gene_symbol),
    n_phenotypes_strict = n_distinct(phenocode),
    n_categories_strict = n_distinct(category),
    genes_strict        = paste(sort(unique(gene_symbol)), collapse = "; "),
    .groups             = "drop"
  )

# 4) Połącz w jedną tabelę
combined_summary <- summary_broad %>%
  full_join(summary_strict, by = c("signature", "random_fdr")) %>%
  arrange(signature, factor(random_fdr, levels = c("0.05","0.1","0.2","ns")))

# 5) Wyświetl wynik
combined_summary %>% as.data.frame() %>% 
  filter(random_fdr != "ns") %>% 
  select(!c(n_categories_broad, n_categories_strict))


# ##############################################################################
# ---- summary for pLoF ----
# ##############################################################################

# 1) Połącz wszystkie sygnatury (p < 0.05) i nadaj rozłączne biny FDR
all_assoc_pLoF <- map_df(signatures, function(sig) {
  grSignature_genebass_pLoF_associations_FDRMonteCarlo$original_genebass_association[[sig]] %>%
    filter(pvalue < 0.05) %>%
    select(-c(q1_pvalue, q3_pvalue, mean_pvalue,
              p30_pvalue, p40_pvalue, p60_pvalue, p80_pvalue, p90_pvalue)) %>%
    mutate(
      signature   = sig,
      fdr_bin     = case_when(
        pvalue < p5_pvalue  ~ "0–0.05",
        pvalue < p10_pvalue ~ "0.05–0.1",
        pvalue < p20_pvalue ~ "0.1–0.2",
        TRUE                ~ ">0.2"
      )
    )
})

# 2) Zbuduj tabelę podsumowującą broad vs strict
summary_pLoF <- all_assoc_pLoF %>%
  group_by(signature, fdr_bin) %>%
  summarise(
    n_var_broad   = dplyr::n(),
    n_genes_broad = n_distinct(gene_symbol),
    n_pheno_broad = n_distinct(phenocode),
    n_var_strict   = sum(pvalue < 1e-5),
    n_genes_strict = n_distinct(gene_symbol[pvalue < 1e-5]),
    n_pheno_strict = n_distinct(phenocode[pvalue < 1e-5]),
    genes_strict   = ifelse(
      n_var_strict > 0,
      paste(sort(unique(gene_symbol[pvalue < 1e-5])), collapse = "; "),
      NA_character_
    ),
    .groups = "drop"
  ) %>%
  arrange(signature,
          factor(fdr_bin, levels = c("0–0.05","0.05–0.1","0.1–0.2",">0.2")))

# 3) Wyświetl wynik
summary_pLoF %>% 
  filter(fdr_bin != ">0.2")
