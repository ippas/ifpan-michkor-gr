GWASCatalogDisGeNETgenebass_GrTissueSystemic_overlapChi2$processed$original_data$df %>% 
  mutate(
    log10_pvalue = -log10(p_value),
    log10_geneOverlap = log10(gene_overlap_count + 1),
    combine_score = log10_pvalue * log10_geneOverlap
  ) %>%
  filter(!(Var1 %in% c(
    "DisGeNET_NA", "genebass_NA", "GWASCatalog_NA",
    "genebass_F9x_Date_F99_first_reported_(mental_disorder,_not_otherwise_specified)"
  ))) %>%
  rowwise() %>%
  mutate(
    matched = list(str_extract_all(Var1, paste(patterns, collapse="|"))[[1]]),
    n_unique_patterns = n_distinct(matched)
  ) %>%
  ungroup() %>%
  filter(n_unique_patterns == 1) %>% 
  mutate(
    icd10_category = map_chr(matched, ~ .x[[1]])
  ) %>%
  # filter(!is.na(regulation), !is.na(icd10_category)) %>%
  filter(!is.na(icd10_category)) %>% 
  select(-c(log2_odds_ratio, fdr_value, fdr)) %>% 
  rename(
    grSignature        = Var2,
    phenotype          = Var1,
    observed_overlap   = gene_overlap_count,
    grSignature_nGenes = Var2_n_genes,
    phenotype_nGenes   = Var1_n_genes,
    mapped_icd10_range = icd10_category
  ) %>% 
  mutate(
    grSignature = case_when(
      str_detect(grSignature, "BloodCellsUp") ~ "bloodUp",
      str_detect(grSignature, "BloodCellsDown") ~ "bloodDown",
      str_detect(grSignature, "LungCellsUp") ~ "lungUp",
      str_detect(grSignature, "LungCellsDown") ~ "lungDown",
      str_detect(grSignature, "NeuralCellsUp") ~ "neuralUp",
      str_detect(grSignature, "NeuralCellsDown") ~ "neuralDown",
      str_detect(grSignature, "global_GR_genes_globalUp") ~ "systemicUp",
      str_detect(grSignature, "global_GR_genes_globalDown") ~ "systemicDown",
      TRUE ~ grSignature
    )
  ) %>% 
  mutate(
    # ---- regulation ----
    regulation = case_when(
      str_detect(grSignature, "Up") ~ "up",
      str_detect(grSignature, "Down") ~ "down",
      TRUE ~ NA_character_
    ),
    
    # ---- tissue ----
    tissue = case_when(
      str_detect(grSignature, "blood") ~ "blood",
      str_detect(grSignature, "lung") ~ "lung",
      str_detect(grSignature, "neural") ~ "neural",
      str_detect(grSignature, "systemic") ~ "systemic",
      TRUE ~ NA_character_
    ),
    
    # ---- signature type ----
    signatureType = case_when(
      tissue == "systemic" ~ "systemic",
      tissue %in% c("blood","lung","neural") ~ "tissue",
      TRUE ~ NA_character_
    )
  ) %>% 
  select(
    mapped_icd10_range,
    phenotype,
    phenotype_nGenes,
    grSignature,
    tissue,
    signatureType,
    regulation,
    grSignature_nGenes,
    expected_overlap,
    observed_overlap,
    odds_ratio,
    chi2,
    p_value,
    combine_score,
    overlap_genes,
    matched,
    n_unique_patterns
  ) %>% 
  select(-matched, -n_unique_patterns) -> preprocessing_summaryTable_grSystemicTissues_icd10F_df

preprocessing_summaryTable_grSystemicTissues_icd10F_df %>% 
  filter(signatureType == "systemic") %>% 
  select(mapped_icd10_range, phenotype, regulation, combine_score) %>% 
  group_by(phenotype) %>% 
  nest() %>% 
  mutate(
    sum_up = map_dbl(data, ~ sum(.x$combine_score[.x$regulation == "up"], na.rm = TRUE)),
    sum_down = map_dbl(data, ~ sum(.x$combine_score[.x$regulation == "down"], na.rm = TRUE)),
    higher_regulation = case_when(
      sum_up > sum_down   ~ "up",
      sum_down > sum_up   ~ "down",
      TRUE                ~ "equal"
    )
  ) %>% 
  unnest() %>% 
  filter(mapped_icd10_range =="F6x") %>% as.data.frame()

preprocessing_summaryTable_grSystemicTissues_icd10F_df %>% 
  filter(signatureType == "systemic") %>% 
  select(mapped_icd10_range, phenotype, regulation, combine_score) %>% 
  group_by(mapped_icd10_range, phenotype, regulation) %>% 
  summarise(sum_cs = sum(combine_score, na.rm = TRUE), .groups = "drop") %>% 
  tidyr::pivot_wider(
    names_from  = regulation,
    values_from = sum_cs,
    values_fill = 0
  ) %>% 
  mutate(
    higher_regulation = case_when(
      up > down  ~ "up",
      down > up  ~ "down",
      TRUE       ~ "equal"
    ),
    score_for_ranking = pmax(up, down)
  ) %>% 
  # filter(higher_regulation %in% c("up", "down")) %>%
  group_by(mapped_icd10_range, higher_regulation) %>% 
  slice_max(order_by = score_for_ranking, n = 3, with_ties = FALSE) %>% 
  mutate(rank = row_number()) %>% 
  ungroup() %>% 
  # filter(mapped_icd10_range == "F9x")
  .$phenotype -> top3grSystemic_phenotypes


top3grSystemic_phenotypes <- c(
  top3grSystemic_phenotypes,
  "DisGeNET_F5x_Pediatric_failure_to_thrive",
  "GWASCatalog_F6x_internet addiction disorder",
  "GWASCatalog_F9x_attention deficit hyperactivity disorder,conduct disorder"
)

# ##############################################################################
# ---- top3 for tissues ----
# ##############################################################################
preprocessing_summaryTable_grSystemicTissues_icd10F_df %>% 
  filter(signatureType == "tissue") %>% 
  select(mapped_icd10_range, phenotype, regulation, combine_score) %>% 
  group_by(mapped_icd10_range, phenotype, regulation) %>% 
  summarise(sum_cs = sum(combine_score, na.rm = TRUE), .groups = "drop") %>% 
  tidyr::pivot_wider(
    names_from  = regulation,
    values_from = sum_cs,
    values_fill = 0
  ) %>% 
  mutate(
    higher_regulation = case_when(
      up > down  ~ "up",
      down > up  ~ "down",
      TRUE       ~ "equal"
    ),
    score_for_ranking = pmax(up, down)
  ) %>% 
  filter(higher_regulation %in% c("up", "down")) %>%
  group_by(mapped_icd10_range, higher_regulation) %>% 
  slice_max(order_by = score_for_ranking, n = 3, with_ties = FALSE) %>% 
  mutate(rank = row_number()) %>% 
  # ungroup() %>% filter(mapped_icd10_range %in% "F7x")
  .$phenotype -> top3grTissues_phenotypes


top3grTissues_phenotypes <- c(
  top3grTissues_phenotypes,
  "DisGeNET_F7x_Coffin-Siris_syndrome"
)


# ##############################################################################
preprocessing_summaryTable_grSystemicTissues_icd10F_df %>% 
  group_by(mapped_icd10_range, signatureType) %>% 
  nest %>%
  mutate(data_signif = map(data, ~ .x %>% filter(p_value < 0.05 & observed_overlap >= 3))) %>% 
    mutate(data_ns     = map(data, ~ .x %>% filter(p_value > 0.05 | observed_overlap < 3))) %>% 
  mutate(signif_cs   = map(data_signif, ~ .x$combine_score)) %>% 
  mutate(ns_cs       = map(data_ns, ~ .x$combine_score)) %>% 
  mutate(n_phenotypes          = map(data, ~ .x$phenotype %>% unique %>% length)) %>% unnest(n_phenotypes) %>% 
  mutate(n_signif_phenotypes   = map(data_signif, ~ .x$phenotype %>% unique %>% length)) %>% unnest(n_signif_phenotypes) %>% 
  mutate(n_signif_associations = map(data_signif, ~ .x %>% nrow)) %>% unnest(n_signif_associations) %>% 
  mutate(mean_signif_cs         = map(data_signif, ~ .x$combine_score %>% mean)) %>% 
  mutate(mean_all_cs            = map(data, ~ .x$combine_score %>% mean)) %>% 
  mutate(sum_signif_cs          = map(data_signif, ~ .x$combine_score %>% sum)) %>% 
  mutate(sum_ns_cs              = map(data_ns, ~ .x$combine_score %>% sum)) %>% 
  unnest(c(mean_signif_cs, mean_all_cs, sum_signif_cs, sum_ns_cs)) %>% 
  mutate(mean_signif_cs         = replace(mean_signif_cs, is.nan(mean_signif_cs), 0)) %>% 
  mutate(data_top3systemic      = map(data, ~ .x %>% filter(phenotype %in% top3grSystemic_phenotypes))) %>% 
  mutate(data_top3Tissues       = map(data, ~ .x %>% filter(phenotype %in% top3grTissues_phenotypes))) %>% 
  mutate(top3Tissues_sumSignif_cs = map(data_top3Tissues, ~.x %>% filter(p_value < 0.05 & observed_overlap >= 3) %>% .$combine_score %>% sum)) %>% 
  mutate(top3Tissues_sumNS_cs   = map(data_top3Tissues, ~ .x %>% filter(p_value > 0.05 | observed_overlap < 3) %>% .$combine_score %>% sum)) %>% 
  mutate(top3Systemic_sumSignif_cs    = map(data_top3systemic, ~.x %>% filter(p_value < 0.05 & observed_overlap >= 3) %>% .$combine_score %>% sum)) %>% 
  mutate(top3Systemic_sumNS_cs  = map(data_top3systemic, ~.x %>% filter(p_value > 0.05 | observed_overlap < 3) %>% .$combine_score %>% sum)) %>% 
  unnest(top3Systemic_sumSignif_cs, top3Systemic_sumNS_cs, top3Tissues_sumSignif_cs, top3Tissues_sumNS_cs) -> preprocessing_summaryTable_grSystemicTissues_icd10F_df
  
preprocessing_summaryTable_grSystemicTissues_icd10F_df %>% 
  mutate(top3Systemic_sumAll_cs = top3Systemic_sumSignif_cs + top3Systemic_sumNS_cs) %>% 
  mutate(sum_all_cs = sum_signif_cs + sum_ns_cs) %>% 
  mutate(nSignif_top3Systemic = map(data_top3systemic, ~ .x %>% filter(p_value < 0.05 & observed_overlap >= 3) %>% nrow)) %>% 
  unnest(nSignif_top3Systemic) %>% 
  mutate(meanSignif_top3Systemic = top3Systemic_sumSignif_cs / nSignif_top3Systemic) %>% 
  mutate(
    meanSignif_top3Systemic = replace(meanSignif_top3Systemic, is.na(meanSignif_top3Systemic), 0)
  ) %>% 
  mutate(meanAll_top3Systemic_cs = case_when(
    signatureType == "tissue" & mapped_icd10_range != "F6x" ~ top3Systemic_sumAll_cs/36,
    signatureType == "tissue" & mapped_icd10_range == "F6x" ~ top3Systemic_sumAll_cs/24,
    signatureType == "systemic" & mapped_icd10_range != "F6x" ~ top3Systemic_sumAll_cs/12,
    signatureType == "systemic" & mapped_icd10_range == "F6x" ~ top3Systemic_sumAll_cs/8
  )) %>% 
  mutate(top3Tissues_sumAll_cs = top3Tissues_sumSignif_cs + top3Tissues_sumNS_cs) %>% 
  mutate(nSignif_top3Tissues = map(data_top3Tissues, ~ .x %>% filter(p_value < 0.05 & observed_overlap >= 3) %>% nrow)) %>% 
  unnest(nSignif_top3Tissues) %>% 
  mutate(meanSignif_top3Tissues = top3Tissues_sumSignif_cs / nSignif_top3Tissues) %>% 
  mutate(
    meanSignif_top3Tissues = replace(meanSignif_top3Tissues, is.na(meanSignif_top3Tissues), 0)
  ) %>% 
  mutate(meanAll_top3Tissues_cs = case_when(
    signatureType == "tissue" & mapped_icd10_range != "F6x" ~ top3Tissues_sumAll_cs/36,
    signatureType == "tissue" & mapped_icd10_range == "F6x" ~ top3Tissues_sumAll_cs/24,
    signatureType == "systemic" & mapped_icd10_range != "F6x" ~ top3Tissues_sumAll_cs/12,
    signatureType == "systemic" & mapped_icd10_range == "F6x" ~ top3Tissues_sumAll_cs/8
  )) -> preprocessing_summaryTable_grSystemicTissues_icd10F_df

preprocessing_summaryTable_grSystemicTissues_icd10F_df %>% 
  select(mapped_icd10_range, n_phenotypes, n_signif_phenotypes)
  
df_plot <- preprocessing_summaryTable_grSystemicTissues_icd10F_df 

df_plot %>% 
  select(-c(data, data_signif, data_ns, signif_cs))

plot_icd10_bar <- function(
    df,
    value_col,
    fill_colors = c(tissue = "#B7966B", systemic = "#552c17"),
    x_range = NULL,
    y_order = NULL,
    zero_one_normalize = FALSE
) {
  
  value_col <- rlang::ensym(value_col)
  
  # opcjonalna normalizacja 0–1
  if (zero_one_normalize) {
    df <- df %>%
      dplyr::group_by(signatureType) %>%
      dplyr::mutate(
        !!value_col := (
          (!!value_col - min(!!value_col, na.rm = TRUE)) /
            (max(!!value_col, na.rm = TRUE) - min(!!value_col, na.rm = TRUE))
        )
      ) %>%
      dplyr::ungroup()
  }
  
  # domyślna kolejność wg systemic
  if (is.null(y_order)) {
    order_levels <- df %>%
      dplyr::filter(signatureType == "systemic") %>%
      dplyr::arrange(!!value_col) %>%
      dplyr::pull(mapped_icd10_range)
  } else {
    order_levels <- y_order
  }
  
  df_plot <- df %>%
    dplyr::mutate(
      mapped_icd10_range = factor(mapped_icd10_range, levels = order_levels),
      signatureType = factor(signatureType, levels = c("tissue", "systemic"))
    )
  
  p <- ggplot2::ggplot(
    df_plot,
    ggplot2::aes(x = !!value_col, y = mapped_icd10_range, fill = signatureType)
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = 0.7),
      width = 0.6,
      color = "black",
      alpha = 0.8
    ) +
    ggplot2::scale_fill_manual(values = fill_colors) +
    ggplot2::labs(
      x = rlang::as_label(value_col),
      y = "ICD-10 category",
      fill = "Signature type"
    ) +
    ggplot2::theme_classic() +
    theme_black_text
  
  if (!is.null(x_range)) {
    p <- p + ggplot2::scale_x_continuous(position = "top", limits = x_range)
  } else {
    p <- p + ggplot2::scale_x_continuous(position = "top")
  }
  
  return(p)
}

plot_icd10_bar(
  preprocessing_summaryTable_grSystemicTissues_icd10F_df,
  mean_signif_cs
) -> p1

plot_icd10_bar(
  preprocessing_summaryTable_grSystemicTissues_icd10F_df,
  mean_all_cs
) -> p2

plot_icd10_bar(
  preprocessing_summaryTable_grSystemicTissues_icd10F_df,
  sum_signif_cs
) -> p3

plot_icd10_bar(
  preprocessing_summaryTable_grSystemicTissues_icd10F_df,
  sum_ns_cs
) -> p4
#
plot_icd10_bar(
  preprocessing_summaryTable_grSystemicTissues_icd10F_df,
  top3Systemic_sumSignif_cs
) -> p5

plot_icd10_bar(
  preprocessing_summaryTable_grSystemicTissues_icd10F_df,
  top3Systemic_sumNS_cs
) -> p6

plot_icd10_bar(
  preprocessing_summaryTable_grSystemicTissues_icd10F_df,
  top3Systemic_sumAll_cs
) -> p7

plot_icd10_bar(
  preprocessing_summaryTable_grSystemicTissues_icd10F_df,
  sum_all_cs
) -> p8

plot_icd10_bar(
  preprocessing_summaryTable_grSystemicTissues_icd10F_df,
  meanAll_top3Systemic_cs
) -> p9

plot_icd10_bar(
  preprocessing_summaryTable_grSystemicTissues_icd10F_df,
  meanSignif_top3Systemic
) -> p10

plot_icd10_bar(
  preprocessing_summaryTable_grSystemicTissues_icd10F_df,
  meanAll_top3Tissues_cs
)  -> p11

plot_icd10_bar(
  preprocessing_summaryTable_grSystemicTissues_icd10F_df,
  meanSignif_top3Tissues
) -> p12


(p1 + p2 + p3 + p4 +
    p5 + p6 + p7 + p8 +
    p9 + p10 + p11 + p12) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")


plot_icd10_bar(
  preprocessing_summaryTable_grSystemicTissues_icd10F_df,
  mean_signif_cs,
  x_range = c(0,5),
  y_order = rev(c("F3x", "F2x", "F0x", "F1x", "F4x", "F8x", "F5x", "F6x", "F7x", "F9x"))
) -> p1

plot_icd10_bar(
  preprocessing_summaryTable_grSystemicTissues_icd10F_df,
  mean_all_cs,
  # x_range = c(0,1),
  y_order = rev(c("F3x", "F2x", "F0x", "F1x", "F4x", "F8x", "F5x", "F6x", "F7x", "F9x"))
) -> p2

plot_icd10_bar(
  preprocessing_summaryTable_grSystemicTissues_icd10F_df,
  mean_all_cs,
  x_range = c(0,5),
  y_order = rev(c("F3x", "F2x", "F0x", "F1x", "F4x", "F8x", "F5x", "F6x", "F7x", "F9x"))
) -> p3


plot_icd10_bar(
  preprocessing_summaryTable_grSystemicTissues_icd10F_df,
  meanAll_top3Systemic_cs,
  x_range = c(0,5),
  y_order = rev(c("F3x", "F2x", "F0x", "F1x", "F4x", "F8x", "F5x", "F6x", "F7x", "F9x"))
) -> p9

plot_icd10_bar(
  preprocessing_summaryTable_grSystemicTissues_icd10F_df,
  meanAll_top3Tissues_cs,
  x_range = c(0,5),
  y_order = rev(c("F3x", "F2x", "F0x", "F1x", "F4x", "F8x", "F5x", "F6x", "F7x", "F9x"))
)  -> p11


svg("/home/mateusz/projects/ifpan-michkor-gr/results_v2/overlap/overlap_grTissues_GWASCatalogDisGeNETgenebass/figures/barplot_grTissuesSystemic_meanCS_v1_04.03.2026.svg",
    width = 8,
    height = 8)

(p1 + p2 + p9 + p11) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

dev.off()
