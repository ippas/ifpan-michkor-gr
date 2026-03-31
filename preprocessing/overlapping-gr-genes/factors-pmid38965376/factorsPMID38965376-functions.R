generate_and_permute_topn_factors_v2 <- function(
    gene_df,
    genes_vector,
    factor_data,             # zawiera kolumnę: signature_name, gene_symbol, pvalue
    n_randomization = 10,
    top_n = 500,
    seed = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
    seeds_vector <- sample(1:1e6, size = n_randomization, replace = FALSE)
  } else {
    seeds_vector <- rep(NA, n_randomization)
  }
  
  original_association <- factor_data %>%
    inner_join(gene_df, by = "gene_symbol") %>%
    select(signature_name, gene_symbol, everything()) %>%
    arrange(pvalue) %>%
    split(.$signature_name) %>%
    lapply(function(df) {
      df %>% head(top_n) %>% mutate(rank = paste0("rank", row_number()))
      
    })
  
  print(head(original_association))
  
  gene_df_list <- set_names(vector("list", n_randomization), paste0("random", seq_len(n_randomization)))
  association_list <- set_names(vector("list", n_randomization), paste0("random", seq_len(n_randomization)))
  
  pb <- txtProgressBar(min = 0, max = n_randomization, style = 3)
  
  for (i in seq_len(n_randomization)) {
    seed_i <- seeds_vector[i]
    if (!is.na(seed_i)) set.seed(seed_i)
    
    random_genes_df <- generate_random_gene_df(
      gene_df = gene_df,
      genes_vector = genes_vector,
      seed = seed_i
    )
    
    assoc <- factor_data %>%
      inner_join(random_genes_df, by = "gene_symbol") %>%
      select(signature_name, gene_symbol, everything()) %>%
      arrange(pvalue) %>%
      split(.$signature_name) %>%
      lapply(function(df) {
        df <- df %>% head(top_n)
        if (nrow(df) > 0) df %>% mutate(rank = paste0("rank", seq_len(dplyr::n())))
        else NULL
      }) %>%
      purrr::compact()
    
    gene_df_list[[i]] <- random_genes_df
    association_list[[i]] <- assoc
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  pvalues_by_signature <- purrr::map(
    names(original_association),
    function(sig_name) {
      purrr::map(association_list, ~ .x[[sig_name]]) %>%
        purrr::keep(~ !is.null(.) && nrow(.) == top_n) %>%
        purrr::map(~ .$pvalue) %>%
        set_names(paste0("random", seq_along(.))) %>%
        as.data.frame() %>%
        mutate(rank = paste0("rank", seq_len(top_n))) %>%
        column_to_rownames(var = "rank")
    }
  ) %>% set_names(names(original_association))
  
  summary_pvalues_by_sig <- purrr::map(
    pvalues_by_signature,
    ~ {
      stats <- apply(.x, 1, function(v) {
        c(
          min = min(v, na.rm = TRUE),
          q1 = quantile(v, 0.25, na.rm = TRUE),
          median = median(v, na.rm = TRUE),
          mean = mean(v, na.rm = TRUE),
          q3 = quantile(v, 0.75, na.rm = TRUE),
          max = max(v, na.rm = TRUE),
          p5 = quantile(v, 0.05, na.rm = TRUE),
          p10 = quantile(v, 0.10, na.rm = TRUE),
          p20 = quantile(v, 0.20, na.rm = TRUE),
          p30 = quantile(v, 0.30, na.rm = TRUE),
          p40 = quantile(v, 0.40, na.rm = TRUE),
          p60 = quantile(v, 0.60, na.rm = TRUE),
          p80 = quantile(v, 0.80, na.rm = TRUE),
          p90 = quantile(v, 0.90, na.rm = TRUE)
        )
      }) %>%
        t() %>%
        as.data.frame() %>%
        tibble::rownames_to_column("rank")
      
      colnames(stats) <- c(
        "rank", "min_pvalue", "q1_pvalue", "median_pvalue", "mean_pvalue", "q3_pvalue", "max_pvalue",
        "p5_pvalue", "p10_pvalue", "p20_pvalue", "p30_pvalue", "p40_pvalue", "p60_pvalue", "p80_pvalue", "p90_pvalue"
      )
      
      stats
    }
  )
  
  original_association <- purrr::imap(
    original_association,
    ~ .x %>% left_join(summary_pvalues_by_sig[[.y]], by = "rank")
  )
  
  list(
    original_association     = original_association,
    gene_df_list             = gene_df_list,
    association_list         = association_list,
    pvalues_by_signature     = pvalues_by_signature,
    summary_pvalues_by_sig   = summary_pvalues_by_sig
  )
}

# ##############################################################################
generate_and_permute_topn_factors_v2 <- function(
    gene_df,
    genes_vector,
    factor_data,             # zawiera kolumnę: signature_name, gene_symbol, pvalue
    n_randomization = 10,
    top_n = 500,
    seed = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
    seeds_vector <- sample(1:1e6, size = n_randomization, replace = FALSE)
  } else {
    seeds_vector <- rep(NA, n_randomization)
  }
  
  original_association <- factor_data %>%
    inner_join(gene_df, by = "gene_symbol") %>%
    select(signature_name, gene_symbol, everything()) %>%
    arrange(pvalue) %>%
    split(.$signature_name) %>%
    lapply(function(df) {
      df %>% head(top_n) %>% mutate(rank = paste0("rank", row_number()))
    })
  
  gene_df_list <- set_names(vector("list", n_randomization), paste0("random", seq_len(n_randomization)))
  association_list <- set_names(vector("list", n_randomization), paste0("random", seq_len(n_randomization)))
  
  pb <- txtProgressBar(min = 0, max = n_randomization, style = 3)
  
  for (i in seq_len(n_randomization)) {
    seed_i <- seeds_vector[i]
    if (!is.na(seed_i)) set.seed(seed_i)
    
    random_genes_df <- generate_random_gene_df(
      gene_df = gene_df,
      genes_vector = genes_vector,
      seed = seed_i
    )
    
    assoc <- factor_data %>%
      inner_join(random_genes_df, by = "gene_symbol") %>%
      select(signature_name, gene_symbol, everything()) %>%
      arrange(pvalue) %>%
      split(.$signature_name) %>%
      lapply(function(df) {
        df <- df %>% head(top_n)
        if (nrow(df) > 0) df %>% mutate(rank = paste0("rank", seq_len(dplyr::n())))
        else NULL
      }) %>%
      purrr::compact()
    
    gene_df_list[[i]] <- random_genes_df
    association_list[[i]] <- assoc
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  pvalues_by_signature <- purrr::map(
    names(original_association),
    function(sig_name) {
      purrr::map(association_list, ~ .x[[sig_name]]) %>%
        purrr::keep(~ !is.null(.) && nrow(.) == top_n) %>%
        purrr::map(~ .$pvalue) %>%
        set_names(paste0("random", seq_along(.))) %>%
        as.data.frame() %>%
        mutate(rank = paste0("rank", seq_len(top_n))) %>%
        column_to_rownames(var = "rank")
    }
  ) %>% set_names(names(original_association))
  
  summary_pvalues_by_sig <- purrr::map(
    pvalues_by_signature,
    ~ {
      stats <- apply(.x, 1, function(v) {
        c(
          min = min(v, na.rm = TRUE),
          q1 = quantile(v, 0.25, na.rm = TRUE),
          median = median(v, na.rm = TRUE),
          mean = mean(v, na.rm = TRUE),
          q3 = quantile(v, 0.75, na.rm = TRUE),
          max = max(v, na.rm = TRUE),
          p5 = quantile(v, 0.05, na.rm = TRUE),
          p10 = quantile(v, 0.10, na.rm = TRUE),
          p20 = quantile(v, 0.20, na.rm = TRUE),
          p30 = quantile(v, 0.30, na.rm = TRUE),
          p40 = quantile(v, 0.40, na.rm = TRUE),
          p60 = quantile(v, 0.60, na.rm = TRUE),
          p80 = quantile(v, 0.80, na.rm = TRUE),
          p90 = quantile(v, 0.90, na.rm = TRUE)
        )
      }) %>%
        t() %>%
        as.data.frame() %>%
        tibble::rownames_to_column("rank")
      
      colnames(stats) <- c(
        "rank", "min_pvalue", "q1_pvalue", "median_pvalue", "mean_pvalue", "q3_pvalue", "max_pvalue",
        "p5_pvalue", "p10_pvalue", "p20_pvalue", "p30_pvalue", "p40_pvalue", "p60_pvalue", "p80_pvalue", "p90_pvalue"
      )
      
      stats
    }
  )
  
  original_association <- purrr::imap(
    original_association,
    ~ .x %>%
      left_join(summary_pvalues_by_sig[[.y]], by = "rank") %>%
      mutate(
        fdr = case_when(
          p5_pvalue      > pvalue ~ "< 0.05",
          p10_pvalue     > pvalue ~ "< 0.1",
          p20_pvalue     > pvalue ~ "< 0.2",
          q1_pvalue      > pvalue ~ "< 0.25",
          p30_pvalue     > pvalue ~ "< 0.3",
          p40_pvalue     > pvalue ~ "< 0.4",
          median_pvalue  > pvalue ~ "< 0.5",
          p60_pvalue     > pvalue ~ "< 0.6",
          q3_pvalue      > pvalue ~ "< 0.75",
          p80_pvalue     > pvalue ~ "< 0.8",
          p90_pvalue     > pvalue ~ "< 0.9",
          TRUE                      ~ "<= 1"
        )
      )
  )
  
  list(
    original_association     = original_association,
    gene_df_list             = gene_df_list,
    association_list         = association_list,
    pvalues_by_signature     = pvalues_by_signature,
    summary_pvalues_by_sig   = summary_pvalues_by_sig
  )
}

process_overlap_heatmap_plot_only <- function(
    factors_df,
    type = c("locus", "geneCenter"),
    window_kb = NULL,
    pvalue_threshold = 1e-5,
    drawing_overlap_threshold = 3,
    pvalue_overlap_threshold = 0.05,
    plot_title_prefix = "",
    col_only_n_genes = TRUE,
    color_scale_range = c(0, 4),
    return_df = FALSE   # << NOWY argument
) {
  type <- match.arg(type)
  
  factor_gene_list <- factors_df %>%
    dplyr::filter(pvalue < pvalue_threshold) %>%
    dplyr::group_by(factor_name) %>%
    dplyr::summarise(genes = list(unique(gene_symbol)), .groups = "drop") %>%
    tibble::deframe()
  
  gene_lists_all <- c(
    split(lite_grSignatures$hgnc_symbol, lite_grSignatures$signature_name),
    factor_gene_list
  )
  
  chi2_results <- perform_chi2_tests(gene_lists_all, total_genes = hgnc_symbols_vector_v110)
  
  rows_to_filter <- unique(lite_grSignatures$signature_name)
  cols_to_filter <- factors_df %>%
    dplyr::filter(pvalue < pvalue_threshold) %>%
    dplyr::pull(factor_name) %>%
    unique()
  
  overlap <- processing_overlap_results(
    data = chi2_results,
    rows_to_filter = rows_to_filter,
    cols_to_filter = cols_to_filter,
    overlap_threshold = 0,
    fdr_threshold = 1,
    genes_list = gene_lists_all
  )
  
  cat("\n=== original_data$df (p_value < ", pvalue_overlap_threshold, ") ===\n", sep = "")
  print(overlap$original_data$df %>% dplyr::filter(p_value < pvalue_overlap_threshold))
  cat("\n")
  
  title_str <- sprintf("%s     %s +/- %s", plot_title_prefix, type,
                       ifelse(!is.null(window_kb), paste0(window_kb, "kb"), "?"))
  
  p <- draw_custom_heatmap_ggplot_v3(
    data_list = overlap,
    title = title_str,
    title_size = 22,
    data_type = "original_data",
    p_thresholds = c(0.05, 0.01),
    color_rects = c("#4C8D05", "#66023C"),
    overlap_threshold = drawing_overlap_threshold,
    apply_filling = FALSE,
    col_only_n_genes = col_only_n_genes,
    col_significant = TRUE,
    color_scale_range = color_scale_range,
    palette = c("white", "#f8dedd", "#f1bcbb", "#edacab", "#e68a89"),
    text_color = "black",
    axis_text_angle = 0,
    axis_text_hjust = 0.5,
    text_size_tile = 4,
    text_size_axis = 16,
    text_size_legend = 18
  )
  
  if (return_df) {
    return(list(plot = p, df = overlap$original_data$df))
  } else {
    return(p)
  }
}


process_overlap_heatmap_plot_only_fisher <- function(factors_df,
                                                     type = c("locus", "geneCenter"),
                                                     window_kb = NULL,
                                                     pvalue_threshold = 1e-5,
                                                     drawing_overlap_threshold = 3,
                                                     fdr_output_threshold = 0.2,
                                                     plot_title_prefix = "",
                                                     col_only_n_genes = TRUE,
                                                     color_scale_range = c(0, 4)) {
  type <- match.arg(type)
  
  factor_gene_list <- factors_df %>%
    dplyr::filter(pvalue < pvalue_threshold) %>%
    dplyr::group_by(factor_name) %>%
    dplyr::summarise(genes = list(unique(gene_symbol)), .groups = "drop") %>%
    tibble::deframe()
  
  gene_lists_all <- c(
    split(lite_grSignatures$hgnc_symbol, lite_grSignatures$signature_name),
    factor_gene_list
  )
  
  fisher_results <- perform_fisher_tests(gene_lists_all, total_genes = hgnc_symbols_vector_v110)
  
  rows_to_filter <- unique(lite_grSignatures$signature_name)
  cols_to_filter <- factors_df %>%
    dplyr::filter(pvalue < pvalue_threshold) %>%
    dplyr::pull(factor_name) %>%
    unique()
  
  overlap <- processing_overlap_results_fisher(
    data = fisher_results,
    rows_to_filter = rows_to_filter,
    cols_to_filter = cols_to_filter,
    overlap_threshold = 0,
    fdr_threshold = 1,
    genes_list = gene_lists_all
  )
  
  cat("\n=== original_data$df (fdr < 0.2) ===\n")
  print(overlap$original_data$df %>% dplyr::filter(fdr < fdr_output_threshold))
  
  title_str <- sprintf("%s     %s +/- %s", plot_title_prefix, type,
                       ifelse(!is.null(window_kb), paste0(window_kb, "kb"), "?"))
  
  p <- draw_custom_heatmap_fisher_ggplot_v3(
    data_list = overlap,
    title = title_str,
    title_size = 22,
    data_type = "original_data",
    p_thresholds = c(0.05, 0.01),
    color_rects = c("#4C8D05", "#66023C"),
    overlap_threshold = drawing_overlap_threshold,
    apply_filling = FALSE,
    col_only_n_genes = col_only_n_genes,
    col_significant = TRUE,
    color_scale_range = color_scale_range,
    palette = c("white", "#f8dedd", "#f1bcbb", "#edacab", "#e68a89"),
    text_color = "black",
    axis_text_angle = 0,
    axis_text_hjust = 0.5,
    text_size_tile = 4,
    text_size_axis = 16,
    text_size_legend = 18
  )
  
  return(p)
}

print_significant_rows <- function(factors_df,
                                                     pvalue_threshold = 1e-5,
                                                     drawing_overlap_threshold = 3,
                                                     pvalue_threshold_fisher = 0.05) {
  
  factor_gene_list <- factors_df %>%
    dplyr::filter(pvalue < pvalue_threshold) %>%
    dplyr::group_by(factor_name) %>%
    dplyr::summarise(genes = list(unique(gene_symbol)), .groups = "drop") %>%
    tibble::deframe()
  
  gene_lists_all <- c(
    split(lite_grSignatures$hgnc_symbol, lite_grSignatures$signature_name),
    factor_gene_list
  )
  
  fisher_results <- perform_fisher_tests(gene_lists_all, total_genes = hgnc_symbols_vector_v110)
  
  rows_to_filter <- unique(lite_grSignatures$signature_name)
  cols_to_filter <- factors_df %>%
    dplyr::filter(pvalue < pvalue_threshold) %>%
    dplyr::pull(factor_name) %>%
    unique()
  
  overlap <- processing_overlap_results_fisher(
    data = fisher_results,
    rows_to_filter = rows_to_filter,
    cols_to_filter = cols_to_filter,
    overlap_threshold = 0,
    fdr_threshold = 1,
    genes_list = gene_lists_all
  )
  
  cat("\n=== original_data$df (fdr < 0.2) ===\n")
  print(overlap$original_data$df %>% dplyr::filter(p_value < pvalue_threshold_fisher) %>% 
          filter(gene_overlap_count >=  drawing_overlap_threshold))
}
