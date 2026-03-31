# ---- old scripts, results from this scipts have bugs ----

# ##############################################################################
# ---- prepare functions ----
# ##############################################################################


generate_random_gene_df <- function(gene_df, genes_vector, seed = NULL) {
  #' Generate a data frame with randomly substituted gene names
  #'
  #' This function takes a data frame with two columns (the first one containing gene names
  #' and the second one containing group/list identifiers) and returns a new data frame
  #' where each unique gene is randomly substituted with a gene from a given vector of genes.
  #' The group assignments remain unchanged.
  #'
  #' @param gene_df A data frame with two columns: `gene_df[[1]]` contains gene names, and `gene_df[[2]]` contains group identifiers.
  #' @param genes_vector A character vector of gene names to use for randomly substituting the original gene names.
  #' @param seed Optional integer for the random number generator seed to ensure reproducible results.
  #'
  #' @return A data frame with the same structure as `gene_df`, where the genes in the first column are randomly substituted.
  #' @examples
  #' gene_df <- data.frame(
  #'   gene = c("GeneA", "GeneB", "GeneC", "GeneA", "GeneB"),
  #'   group = c("List1", "List2", "List2", "List1", "List2")
  #' )
  #' genes_vector <- c("GeneX", "GeneY", "GeneZ", "GeneW")
  #' generate_random_gene_df(gene_df, genes_vector, seed = 123)
  #'
  
  # Validate input is a data frame with two columns
  if (!is.data.frame(gene_df)) {
    stop("Input must be a data frame.")
  }
  if (ncol(gene_df) != 2) {
    stop("Input data frame must have exactly 2 columns: one for genes and one for groups.")
  }
  
  # Set the seed if provided to ensure reproducibility
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Extract unique original genes
  original_genes <- unique(gene_df[[1]])
  
  # Sample the same number of genes from the provided vector
  random_genes <- sample(genes_vector, length(original_genes))
  
  # Create a lookup vector mapping original genes to new random genes
  gene_replacement <- setNames(random_genes, original_genes)
  
  # Replace genes in the data frame using the mapping
  gene_df[[1]] <- gene_replacement[gene_df[[1]]]
  
  # Return the updated data frame
  return(gene_df)
}



generate_and_permutate_topn_v2 <- function(
    genebass_data,
    genes_df,
    top_n = 1000,
    n_permutation = 100,
    seed = NULL,
    return_full_results = FALSE,
    return_median_pvalue = TRUE,
    return_mapper_genes = FALSE
) {
  # Ustawienie seed i generacja wektora seedów dla permutacji
  if (!is.null(seed)) {
    set.seed(seed)
    seed_vector <- sample.int(.Machine$integer.max, n_permutation)
  } else {
    seed_vector <- rep(NA, n_permutation)
  }
  
  signatures <- unique(genes_df$signature_name)
  all_genes <- unique(genes_df$gene_symbol)
  
  pvalue_matrices <- setNames(
    vector("list", length(signatures)),
    signatures
  )
  
  for (sig in signatures) {
    pb <- txtProgressBar(min = 0, max = n_permutation, style = 3)
    perm_pvalues <- matrix(NA_real_, nrow = top_n, ncol = n_permutation)
    full_results_list <- list()
    mapper_list <- list()  # przechowujemy mapowania jeśli włączone
    
    for (i in seq_len(n_permutation)) {
      if (!is.na(seed_vector[i])) set.seed(seed_vector[i])
      
      # Stwórz globalne mapowanie dla tej permutacji
      mapping <- sample(all_genes, length(all_genes))
      names(mapping) <- all_genes
      
      # Jeśli włączone, zapisujemy mapowanie dla tej permutacji
      if (return_mapper_genes) {
        mapper_list[[i]] <- mapping
      }
      
      # Podmień geny w całym genes_df
      permuted_df <- genes_df %>%
        dplyr::mutate(
          gene_symbol = mapping[as.character(gene_symbol)]
        )
      
      # Wybierz podzbiór dla tej sygnatury
      sig_df <- permuted_df[permuted_df$signature_name == sig, ]
      merged_df <- merge(
        genebass_data,
        sig_df,
        by = "gene_symbol"
      )
      top_df <- merged_df[order(merged_df$pvalue), ]
      top_vals <- head(top_df$pvalue, top_n)
      perm_pvalues[, i] <- top_vals
      
      if (return_full_results) {
        full_results_list[[i]] <- top_df
      }
      setTxtProgressBar(pb, i)
    }
    close(pb)
    
    rownames(perm_pvalues) <- paste0("rank", seq_len(top_n))
    colnames(perm_pvalues) <- paste0("random", seq_len(n_permutation))
    
    median_df <- NULL
    if (return_median_pvalue) {
      median_df <- data.frame(
        rank = rownames(perm_pvalues),
        median_pvalue = apply(perm_pvalues, 1, median, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
    
    result_list <- list(
      pvalues_matrix = perm_pvalues
    )
    if (!is.null(median_df)) {
      result_list$median_df <- median_df
    }
    if (return_full_results) {
      names(full_results_list) <- paste0("random", seq_len(n_permutation))
      result_list$full_results <- full_results_list
    }
    if (return_mapper_genes) {
      names(mapper_list) <- paste0("random", seq_len(n_permutation))
      result_list$mapper_genes <- mapper_list
    }
    pvalue_matrices[[sig]] <- result_list
  }
  
  return(pvalue_matrices)
}

make_pvalue_hist <- function(
    pvalues, 
    breaks = 10,           # domyślnie 10
    main = NULL, 
    col = "gray"
) {
  # transformacja
  logp <- -log10(pvalues)
  med <- median(logp, na.rm = TRUE)
  
  # histogram
  hist(
    logp,
    breaks = breaks,
    main = if (is.null(main)) "Histogram -log10(p-values)" else main,
    xlab = "-log10(p)",
    col = col,
    border = "black"
  )
  
  # linia dla mediany
  abline(v = med, col = "red", lwd = 2)
  
  # tekst po prawej stronie linii
  text(
    x = med,
    y = par("usr")[4] * 0.9,
    labels = paste0("median = ", round(med, 2)),
    col = "red",
    pos = 4,      # tekst po prawej
    offset = 0.2
  )
}
# ##############################################################################
# ---- prepare data ----
# ##############################################################################


lite_grSignatures <- gr_genes_signatures_multi_approach_df %>% 
  filter(signature_name %in% c("universal_up", "universal_down", "brain_up", "brain_down")) %>% 
  dplyr::select(c(hgnc_symbol, signature_name)) %>% 
  set_colnames(c("gene_symbol", "signature_name")) %>% 
  mutate(signature_name = case_when(
    signature_name == "universal_up"   ~ "metasignature_up",
    signature_name == "universal_down" ~ "metasignature_down",
    TRUE ~ signature_name
  ))

genebass_mentalHealth_skat <-read.delim("data/genebass/mentalHealth_Pvalue_SKAT_0.05.tsv.bgz")




# ##############################################################################
# ---- permutation test ----
# ##############################################################################
pvalue_matrices_v2 <- generate_and_permutate_topn_v2(
  genebass_data = genebass_mentalHealth_skat,
  genes_df = lite_grSignatures,
  top_n = 500,
  n_permutation = 100,
  seed = 123,
  return_full_results = TRUE,
  return_median_pvalue = TRUE,
  return_mapper_genes = TRUE
)

pvalue_matrices_v2$metasignature_down$pvalues_matrix %>% .[1, ] %>% log10 %>% {. * (-1)} %>% hist
pvalue_matrices_v2$metasignature_down$pvalues_matrix %>% .[1, ] %>% log10 %>% {. * (-1)} %>% length()
pvalue_matrices_v2$metasignature_down$pvalues_matrix %>% .[1, ]  %>% sort

pvalue_matrices_v2$metasignature_down$full_results$random1 %>% head(1)
pvalue_matrices_v2$metasignature_down$full_results$random9 %>% head(1)
pvalue_matrices_v2$metasignature_down$full_results$random14 %>% head(1)


pvalue_matrices_v2$metasignature_down$mapper_genes$random1 %>% as.character() %>% sort
pvalue_matrices_v2$metasignature_down$mapper_genes$random14 %>% as.character() %>% sort
pvalue_matrices_v2$metasignature_down$mapper_genes$random100 %>% as.character() %>% sort

pvalue_matrices_v2$metasignature_down$mapper_genes$random99 %>% as.character() %>% sort
pvalue_matrices_v2$metasignature_up$mapper_genes %>% unlist %>% unname() %>% unique



pvalue_matrices_v2$brain_down$median_df 

pvalue_matrices_v2$metasignature_up$mapper_genes$random1["LPL"]

pvalue_matrices_v2$brain_down$mapper_genes$random1["LPL"]


grSignaturesLite_association_genebassSKATp05 %>% 
  lapply(., function(x){
    x %>% arrange(pvalue) %>% 
      head(500) %>% 
      mutate(rank_pvalue = paste0("rank", seq(1:500))) %>% 
      select(rank_pvalue, everything())
  }) -> grSignaturesLite_association_genebassSKATtop500


grSignaturesLite_association_genebassSKATtop500 <- mapply(
  function(obs_df, perm_list) {
    # Wyciągamy median_df z perm_list i zmieniamy kolumnę
    median_df <- perm_list$median_df %>%
      dplyr::rename(median_perm_pvalue = median_pvalue) %>%
      dplyr::mutate(rank_pvalue = as.character(rank)) %>%
      dplyr::select(rank_pvalue, median_perm_pvalue)
    
    obs_df %>%
      dplyr::left_join(median_df, by = "rank_pvalue")
  },
  grSignaturesLite_association_genebassSKATtop500,
  pvalue_matrices_v2,
  SIMPLIFY = FALSE
)

grSignaturesLite_association_genebassSKATtop500$metasignature_down %>% 
  filter(pvalue == median_perm_pvalue)


grSignaturesLite_association_genebassSKATtop500 <- lapply(
  grSignaturesLite_association_genebassSKATtop500,
  function(df) {
    df %>%
      dplyr::mutate(
        p_less_than_permuted = pvalue < median_perm_pvalue
      ) %>%
      dplyr::select(
        rank_pvalue, pvalue, median_perm_pvalue, p_less_than_permuted,
        dplyr::everything()
      ) %>%
      dplyr::select(
        signature_name, gene_id, gene_symbol, in_GWASCatalog,
        annotation, phenocode, description, category, pvalue,
        median_perm_pvalue, p_less_than_permuted, beta, pvalue_test
      )
  }
)


pvalue_matrices_v2$metasignature_up$pvalues_matrix[1:50, ] %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "rank") %>% 
  pivot_longer(
    cols = -rank,
    names_to = "permutation",
    values_to = "p"
  ) %>%
  mutate(
    # wyciągamy numer z "rank1" -> 1
    rank_num = as.numeric(str_remove(rank, "rank")),
    log10p = -log10(p)
  ) %>%
  ggplot(aes(x = factor(rank_num), y = log10p)) +
  geom_boxplot(outlier.shape = NA) +
  labs(
    title = "-log10(p) dla permutowanych list (50 rang)",
    x = "Ranga",
    y = "-log10(p)"
  ) +
  scale_y_continuous(
    breaks = seq(0, 15, 1) # linie co 1 w przedziale 0-15
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid.major.x = element_blank(),  # pionowych nie chcemy
    panel.grid.minor = element_blank()      # minor niepotrzebne
  )


df_plot <- pvalue_matrices_v2$metasignature_up$pvalues_matrix[1:50, ] %>%
  as.data.frame() %>%
  rownames_to_column(var = "rank") %>%
  pivot_longer(
    cols = -rank,
    names_to = "permutation",
    values_to = "p"
  ) %>%
  mutate(
    rank_num = as.numeric(str_remove(rank, "rank")),
    log10p = -log10(p)
  )

# Podsumowanie per ranga
df_summary <- df_plot %>%
  group_by(rank_num) %>%
  summarise(
    min = min(log10p),
    q1 = quantile(log10p, 0.25),
    q3 = quantile(log10p, 0.75),
    max = max(log10p),
    range = max - min,
    iqr = q3 - q1
  ) %>%
  ungroup()

ggplot(df_plot, aes(x = factor(rank_num), y = log10p)) +
  geom_boxplot(outlier.shape = NA) +
  geom_text(
    data = df_summary,
    aes(
      x = factor(rank_num),
      y = max + 0.1,
      label = paste0(
        "min=", round(min,1),
        ", max=", round(max,1),
        ", Δ=", round(range,1),
        ", IQR=", round(iqr,1)
      )
    ),
    size = 2.5,
    angle = 90,
    hjust = 0
  ) +
  scale_y_continuous(breaks = seq(0, 15, 1)) +
  labs(
    title = "-log10(p) dla permutowanych list (50 rang)",
    x = "Ranga",
    y = "-log10(p)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


grSignaturesLite_association_genebassSKATtop500$metasignature_down %>% 
  mutate(rank_pvalue = paste0("rank", seq(1:500))) %>% 
  filter(pvalue == median_perm_pvalue) %>% .$rank_pvalue

grSignaturesLite_association_genebassSKATtop500 %>% 
  lapply(., function(x){
    x %>% 
      mutate(log_pvalue = (-1)*log10(pvalue)) %>% 
      mutate(log_median_pvalue = (-1)*log10(median_perm_pvalue)) %>% 
      mutate(diff_pvalue = log_pvalue - log_median_pvalue)
  }) %>% 
  .$brain_up %>% .$diff_pvalue %>% hist


grSignaturesLite_association_genebassSKATtop500 %>% 
  lapply(function(x) {
    x %>% mutate(
      log_pvalue = -log10(pvalue),
      log_median_pvalue = -log10(median_perm_pvalue),
      diff_pvalue = abs(log_pvalue - log_median_pvalue)
    )
  }) %>% 
  bind_rows(.id = "category") %>% 
  ggplot(aes(x = diff_pvalue)) +
  geom_histogram(bins = 50, color = "black", fill = "lightgray") +
  facet_wrap(~ category, ncol = 2) +
  theme_minimal() +
  labs(
    x = "|log10(p) - log10(median perm p)|",
    y = "Count",
    title = "Rozkłady różnic dla poszczególnych sygnatur, \n23.06.2025r, stare wyniki"
  )


grSignaturesLite_association_genebassSKATtop500 %>% 
  lapply(function(x) {
    x %>% mutate(
      log_pvalue = -log10(pvalue),
      log_median_pvalue = -log10(median_perm_pvalue),
      diff_pvalue = abs(log_pvalue - log_median_pvalue)
    )
  }) %>% 
  bind_rows(.id = "category") %>% 
  ggplot(aes(x = diff_pvalue)) +
  geom_histogram(bins = 100, color = "black", fill = "lightgray") +
  facet_wrap(~ category, ncol = 2) +
  coord_cartesian(xlim = c(0, 1)) +   # <- max do 2
  theme_minimal() +
  labs(
    x = "log10(p) - log10(median perm p)",
    y = "Count",
    title = "Rozkłady różnic dla poszczególnych sygnatur, \n23.06.2025r, stare wyniki, oś x od 0 do 1"
  )

grSignaturesLite_association_genebassSKATtop500 %>% 
  lapply(function(x) {
    x %>% mutate(
      log_pvalue = -log10(pvalue),
      log_median_pvalue = -log10(median_perm_pvalue),
      diff_pvalue = log_pvalue - log_median_pvalue
    )
  }) %>% 
  sapply(function(df) e1071::kurtosis(df$diff_pvalue))

# 
# # ##############################################################################
# # ---- save to xlsx ----
# # ##############################################################################
# wb <- createWorkbook()
# 
# for (sheetname in c("metasignature_up", "metasignature_down", "brain_up", "brain_down")) {
#   df <- grSignaturesLite_association_genebassSKATtop500[[sheetname]]
#   
#   addWorksheet(wb, sheetName = sheetname)
#   writeData(wb, sheet = sheetname, x = df, withFilter = TRUE)
#   freezePane(wb, sheet = sheetname, firstRow = TRUE)
# }
# 
# # 5. Zapisz workbook
# saveWorkbook(wb, file = "data/genebass/grSignaturesLite_association_genebassSKATtop500.xlsx", overwrite = TRUE)
