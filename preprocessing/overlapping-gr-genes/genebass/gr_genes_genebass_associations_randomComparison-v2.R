# ##############################################################################
# ---- validation of permutations ----
# ##############################################################################

# ##############################################################################
# ---- read data ----
# ##############################################################################
genebass_mentalHealth_skat_all <-read.delim("data/genebass/genebass_mental_health_all_SKAT.tsv.bgz")


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


plot_permutation_pvalues_by_rank_v3 <- function(
    signature_data,
    ranks_to_plot = NULL,
    x_axis_limit = 12,
    title_text_size = 16,
    subtitle_text_size = 12,
    axis_title_size = 14,
    axis_text_x_size = 12,
    axis_text_y_size = 12,
    legend_text_size = 10,
    legend_title_size = 12,
    legend_title = "P-value summary",
    width_rank = 3,
    width_gene = 8,
    width_annot = 12,
    width_desc = 40,
    use_mono_font = TRUE,
    pvalue_types_to_plot = c("pvalue", "median_pvalue", "mean_pvalue"),
    point_size = 2,
    highlight_point_size = 3,
    highlight_point_shape = 8,
    pvalue_color_map_override = NULL
) {
  # Filter by requested ranks if provided
  if (!is.null(ranks_to_plot)) {
    if (!is.character(ranks_to_plot)) {
      stop("ranks_to_plot must be a character vector (e.g. c('rank1','rank2')).")
    }
    my_data_filtered <- signature_data %>%
      dplyr::filter(rank %in% ranks_to_plot)
  } else {
    my_data_filtered <- signature_data
  }
  
  if (nrow(my_data_filtered) == 0) {
    warning("No data left after filtering by ranks_to_plot. Returning NULL.")
    return(NULL)
  }
  
  # Reshape to long format with -log10(pvalue)
  my_data_long <- my_data_filtered %>%
    dplyr::select(c("rank", pvalue_types_to_plot)) %>%
    tidyr::pivot_longer(
      cols = tidyselect::all_of(pvalue_types_to_plot),
      names_to = "pvalue_type",
      values_to = "pvalue_raw"
    ) %>%
    dplyr::mutate(pvalue_value = -log10(pmax(pvalue_raw, .Machine$double.xmin)))
  
  # Build rank labels (mono or plain)
  rank_labels_df <- my_data_filtered %>%
    dplyr::select(rank, gene_symbol, annotation, description) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      numeric_rank = as.numeric(gsub("rank", "", rank)),
      rank_label = if (use_mono_font) {
        sprintf(paste0("%-", width_rank, "s %-", width_gene, "s %-", width_annot, "s %-", width_desc, "s"),
                numeric_rank, gene_symbol, annotation, description)
      } else {
        paste(numeric_rank, gene_symbol, annotation, description, sep = " --- ")
      },
      rank_label = factor(rank_label, levels = rank_label)
    )
  
  my_data_long <- my_data_long %>%
    dplyr::left_join(rank_labels_df, by = "rank")
  
  # Default and overridden color maps
  default_color_map <- c(
    "pvalue"        = "blue",
    "min_pvalue"    = "#E41A1C",
    "q1_pvalue"     = "#FF7F00",
    "median_pvalue" = "#984EA3",
    "mean_pvalue"   = "#377EB8",
    "q3_pvalue"     = "#4DAF4A",
    "max_pvalue"    = "#A65628",
    "p10_pvalue"    = "#B3CDE3",
    "p20_pvalue"    = "#CCEBC5",
    "p30_pvalue"    = "#DECBE4",
    "p40_pvalue"    = "#FED9A6",
    "p60_pvalue"    = "#FFFFCC",
    "p80_pvalue"    = "#E5D8BD",
    "p90_pvalue"    = "#FDDAEC"
  )
  pvalue_label_map <- c(
    "pvalue"        = "P (GR-dependent gene)",
    "min_pvalue"    = "Min P (random)",
    "q1_pvalue"     = "25th P percentile (random)",
    "median_pvalue" = "Median P (random)",
    "mean_pvalue"   = "Mean P (random)",
    "q3_pvalue"     = "75th P percentile (random)",
    "max_pvalue"    = "Max P (random)",
    "p5_pvalue"     = "5th percentile P (random)",
    "p10_pvalue"    = "10th percentile P (random)",
    "p20_pvalue"    = "20th percentile P (random)",
    "p30_pvalue"    = "30th percentile P (random)",
    "p40_pvalue"    = "40th percentile P (random)",
    "p60_pvalue"    = "60th percentile P (random)",
    "p80_pvalue"    = "80th percentile P (random)",
    "p90_pvalue"    = "90th percentile P (random)"
  )
  color_map_final <- if (!is.null(pvalue_color_map_override)) {
    pvalue_color_map_override
  } else {
    default_color_map[pvalue_types_to_plot]
  }
  
  # Build the plot
  p <- ggplot2::ggplot(my_data_long, ggplot2::aes(x = pvalue_value, y = rank_label)) +
    ggplot2::geom_point(ggplot2::aes(color = pvalue_type), size = point_size) +
    ggplot2::geom_point(
      data = my_data_long %>% dplyr::filter(pvalue_type == "pvalue"),
      ggplot2::aes(x = pvalue_value, y = rank_label),
      color = color_map_final[["pvalue"]],
      size = highlight_point_size,
      shape = highlight_point_shape
    ) +
    ggplot2::xlim(0, x_axis_limit) +
    ggplot2::labs(
      title    = paste0("Distribution of P-values for signature: ", unique(signature_data$signature_name)),
      subtitle = "With summary statistics derived from random sampling and P-value for the GR-dependent gene",
      x        = "-log10(P-value)",
      y        = NULL,
      color    = legend_title
    ) +
    ggplot2::scale_y_discrete(
      breaks   = rank_labels_df$rank_label,
      limits   = rev(rank_labels_df$rank_label),
      position = "right"
    ) +
    ggplot2::scale_color_manual(
      values = color_map_final,
      labels = pvalue_label_map[pvalue_types_to_plot]
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title           = ggplot2::element_text(hjust = 0.5, face = "bold", size = title_text_size),
      plot.subtitle        = ggplot2::element_text(hjust = 0.5, size = subtitle_text_size),
      
      # CZARNE podpisy i ticki osi
      axis.title.x         = ggplot2::element_text(color = "black", size = axis_title_size, face = "bold"),
      axis.title.y         = ggplot2::element_text(color = "black", size = axis_title_size, face = "bold"),
      axis.text.x          = ggplot2::element_text(color = "black", size = axis_text_x_size),
      axis.text.y.right    = ggplot2::element_text(color = "black", size = axis_text_y_size, family = if (use_mono_font) "mono" else ""),
      axis.ticks.y.right   = ggplot2::element_line(),
      
      panel.grid.major.y   = ggplot2::element_blank(),
      panel.grid.minor.y   = ggplot2::element_blank(),
      legend.position      = "bottom",
      legend.text          = ggplot2::element_text(size = legend_text_size),
      legend.title         = ggplot2::element_text(size = legend_title_size, face = "bold")
    )
  
  return(p)
}




generate_and_permute_topn_v3 <- function(
    gene_df,
    genes_vector,
    genebass_data,
    n_randomization = 10,
    top_n = 500,
    seed = NULL
) {
  
  #' generate_and_permute_topn_v3
  #'
  #' This function extracts top-ranked gene-phenotype associations for a specified list of signatures, 
  #' then generates randomized gene lists for a number of permutations and computes summary statistics of permuted p-values.
  #' Finally, it appends summary statistics to the original top-ranked associations.
  #'
  #' @param gene_df Data frame of genes (must include `gene_id` and `gene_symbol`)
  #' @param genes_vector Character vector of all possible gene symbols
  #' @param genebass_data Data frame with Genebass summary statistics (must include `pvalue`, `signature_name`, and `gene_symbol`)
  #' @param n_randomization Number of randomizations to perform
  #' @param top_n Number of top hits to extract per signature
  #' @param seed Optional seed for reproducible random sampling
  #'
  #' @return A list with:
  #' \itemize{
  #'   \item original_genebass_association — named list of top_n hits per signature in the original data
  #'   \item gene_df_list — list of randomly generated gene data frames
  #'   \item genebass_association_list — list of Genebass associations per randomized gene list
  #'   \item pvalues_by_signature — list of p-value matrices per signature across all randomizations
  #'   \item summary_pvalues_by_sig — list of summary statistics for permuted p-values per signature
  #' }
  
  if (!is.null(seed)) {
    set.seed(seed)
    seeds_vector <- sample(1:1e6, size = n_randomization, replace = FALSE)
  } else {
    seeds_vector <- rep(NA, n_randomization)
  }
  
  original_genebass_association <- genebass_data %>%
    inner_join(gene_df, by = "gene_symbol") %>%
    select(signature_name, gene_id, gene_symbol, everything()) %>%
    arrange(pvalue) %>%
    select(-pvalue_threshold, -heritability, -pvalue_test) %>%
    split(.$signature_name) %>%
    lapply(function(df) {
      df %>% head(top_n) %>% mutate(rank = paste0("rank", seq_len(top_n)))
    })
  
  gene_df_list <- set_names(vector("list", n_randomization), paste0("random", seq_len(n_randomization)))
  genebass_association_list <- set_names(vector("list", n_randomization), paste0("random", seq_len(n_randomization)))
  
  pb <- txtProgressBar(min = 0, max = n_randomization, style = 3)
  
  for (i in seq_len(n_randomization)) {
    name_random <- paste0("random", i)
    seed_i      <- seeds_vector[i]
    if (!is.na(seed_i)) set.seed(seed_i)
    
    random_genes_df <- generate_random_gene_df(
      gene_df      = gene_df,
      genes_vector = genes_vector,
      seed         = if (!is.na(seed_i)) seed_i else NULL
    )
    
    assoc <- genebass_data %>%
      inner_join(random_genes_df, by = "gene_symbol") %>%
      select(signature_name, gene_id, gene_symbol, everything()) %>%
      arrange(pvalue) %>%
      select(-pvalue_threshold, -heritability, -pvalue_test) %>%
      split(.$signature_name) %>%
      lapply(function(df) {
        df %>% head(top_n) %>% mutate(rank = paste0("rank", seq_len(top_n)))
      })
    
    gene_df_list[[name_random]] <- random_genes_df
    genebass_association_list[[name_random]] <- assoc
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  pvalues_by_signature <- purrr::map(
    names(original_genebass_association),
    function(sig_name) {
      purrr::map(genebass_association_list, ~ .x[[sig_name]]$pvalue) %>%
        set_names(names(genebass_association_list)) %>%
        as.data.frame() %>%
        mutate(rank = paste0("rank", seq_len(top_n))) %>%
        column_to_rownames(var = "rank")
    }
  ) %>% set_names(names(original_genebass_association))
  
  summary_pvalues_by_sig <- purrr::map(
    pvalues_by_signature,
    ~ {
      basic_stats <- apply(.x, 1, function(v) {
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
        tibble::rownames_to_column(var = "rank")
      
      colnames(basic_stats) <- c(
        "rank",
        "min_pvalue", "q1_pvalue", "median_pvalue", "mean_pvalue", "q3_pvalue", "max_pvalue",
        "p5_pvalue", "p10_pvalue", "p20_pvalue", "p30_pvalue", "p40_pvalue", "p60_pvalue", "p80_pvalue", "p90_pvalue"
      )
      
      basic_stats
    }
  )
  
  original_genebass_association <- purrr::imap(
    original_genebass_association,
    ~ .x %>% dplyr::left_join(summary_pvalues_by_sig[[.y]], by = "rank")
  )
  
  list(
    original_genebass_association = original_genebass_association,
    gene_df_list                  = gene_df_list,
    genebass_association_list     = genebass_association_list,
    pvalues_by_signature          = pvalues_by_signature,
    summary_pvalues_by_sig        = summary_pvalues_by_sig
  )
}



# ##############################################################################
# ---- 1. analysis  ----
# ##############################################################################

# ---- all variants ----
generate_and_permute_topn_v3(gene_df = lite_grSignatures %>% 
                               select(-signature_derivation) %>% 
                               set_colnames(c("gene_symbol", "signature_name")),
                             genes_vector = hgnc_symbols_vector_v110,
                             genebass_data = genebass_mentalHealth_skat_all,
                             n_randomization = 1000,
                             top_n = 1000
) -> grSignature_genebass_associations_FDRMonteCarlo

grSignature_genebass_associations_FDRMonteCarlo$original_genebass_association$brain_up %>% 
  mutate(fdr_0.1 = ifelse(pvalue < p10_pvalue, T, F)) %>% 
  filter(fdr_0.1) %>% 
  .$gene_symbol %>% unique %>% length()

grSignature_genebass_associations_FDRMonteCarlo$original_genebass_association %>% 
  lapply(., function(x){
    x %>% 
      mutate(fdr_0.1 = ifelse(pvalue < p10_pvalue, T, F)) %>% 
      filter(fdr_0.1)
  }) %>% 
  bind_rows() %>% 
  filter(signature_name == "brain_up") %>% .$gene_symbol %>% unique %>% cat(sep = "\n")






# ---- pLoF ----
generate_and_permute_topn_v3(gene_df = lite_grSignatures,
                             genes_vector = hgnc_symbols_vector_v110,
                             genebass_data = genebass_mentalHealth_skat_all %>% filter(annotation == "pLoF"),
                             n_randomization = 1000,
                             top_n = 1000
) -> grSignature_genebass_pLoF_associations_FDRMonteCarlo

grSignature_genebass_pLoF_associations_FDRMonteCarlo$original_genebass_association$brain_up %>% 
  mutate(fdr_0.05 = ifelse(pvalue < p5_pvalue, T, F)) %>%
  mutate(fdr_0.1 = ifelse(pvalue < p10_pvalue, T, F)) %>% 
  mutate(fdr_0.2 = ifelse(pvalue < p20_pvalue, T, F)) %>% 
  mutate(fdr_0.3 = ifelse(pvalue < p30_pvalue, T, F)) %>% head
  filter(fdr_0.3) %>% head



# ---- synonymous ----
  generate_and_permute_topn_v3(gene_df = lite_grSignatures,
                               genes_vector = hgnc_symbols_vector_v110,
                               genebass_data = genebass_mentalHealth_skat_all %>% filter(annotation == "synonymous"),
                               n_randomization = 1000,
                               top_n = 1000
  ) -> grSignature_genebass_synonymous_associations_FDRMonteCarlo

  
grSignature_genebass_synonymous_associations_FDRMonteCarlo$original_genebass_association$brain_up %>% 
  mutate(fdr_0.05 = ifelse(pvalue < p5_pvalue, T, F)) %>%
  mutate(fdr_0.1 = ifelse(pvalue < p10_pvalue, T, F)) %>% 
  mutate(fdr_0.2 = ifelse(pvalue < p20_pvalue, T, F)) %>% 
  mutate(fdr_0.3 = ifelse(pvalue < p30_pvalue, T, F)) %>% 
  filter(fdr_0.2)





# ##############################################################################
# ---- visualization ----
# ##############################################################################




genebass_mentalHealth_skat %>% filter(annotation == "pLoF")


randomization_results_plof <- generate_and_permute_topn_v2(gene_df = lite_grSignatures,
                                                      genes_vector = hgnc_symbols_vector_v110,
                                                      genebass_data = genebass_mentalHealth_skat,
                                                      n_randomization = 1000
)

# ##############################################################################
# ---- save to xlsx ----
# ##############################################################################
setNames(
  lapply(c("metasignature_up", "metasignature_down", "brain_up", "brain_down"), function(sheetname) {
    grSignature_genebass_associations_FDRMonteCarlo$original_genebass_association[[sheetname]] %>%
      filter(pvalue < 0.05) %>%
      select(!c(q1_pvalue, q3_pvalue, mean_pvalue, p30_pvalue, p40_pvalue, p60_pvalue, p80_pvalue, p90_pvalue)) %>%
      mutate(
        signature_name = sheetname,
        category_simple = str_replace(category, "Online follow-up > Mental health > ", ""),
        random_fdr = case_when(
          pvalue < p5_pvalue  ~ "0.05",
          pvalue < p10_pvalue ~ "0.1",
          pvalue < p20_pvalue ~ "0.2",
          TRUE                ~ "ns."
        )
      ) %>%
      select(c(signature_name, gene_id, gene_symbol, annotation, phenocode, description, category, category_simple,
               pvalue, rank, min_pvalue, median_pvalue, max_pvalue, p5_pvalue, p10_pvalue, p20_pvalue, random_fdr))
  }),
  c("metasignature_up", "metasignature_down", "brain_up", "brain_down")
) %>% 
  save_gr_signature_associations_to_xlsx(
    output_file = "data/genebass/grSignaturesLite_associtionSKAT05_withRandomFDR-01.07.2025.xlsx"
  )


setNames(
  lapply(c("metasignature_up", "metasignature_down", "brain_up", "brain_down"), function(sheetname) {
    grSignature_genebass_pLoF_associations_FDRMonteCarlo$original_genebass_association[[sheetname]] %>%
      filter(pvalue < 0.05) %>%
      select(!c(q1_pvalue, q3_pvalue, mean_pvalue, p30_pvalue, p40_pvalue, p60_pvalue, p80_pvalue, p90_pvalue)) %>%
      mutate(
        signature_name = sheetname,
        category_simple = str_replace(category, "Online follow-up > Mental health > ", ""),
        random_fdr = case_when(
          pvalue < p5_pvalue  ~ "0.05",
          pvalue < p10_pvalue ~ "0.1",
          pvalue < p20_pvalue ~ "0.2",
          TRUE                ~ "ns."
        )
      ) %>%
      select(c(signature_name, gene_id, gene_symbol, annotation, phenocode, description, category, category_simple,
               pvalue, rank, min_pvalue, median_pvalue, max_pvalue, p5_pvalue, p10_pvalue, p20_pvalue, random_fdr))
  }),
  c("metasignature_up", "metasignature_down", "brain_up", "brain_down")
) %>% 
  save_gr_signature_associations_to_xlsx(
    output_file = "data/genebass/grSignaturesLite_associtionSKAT05pLoF_withRandomFDR-01.07.2025.xlsx"
  )
