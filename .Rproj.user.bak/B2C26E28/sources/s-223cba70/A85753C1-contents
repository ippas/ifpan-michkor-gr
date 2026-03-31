# ##############################################################################
# ---- prepare data ----
# ##############################################################################

lite_grSignatures
hgnc_symbols_vector_v110

factors_rsidGenes_P1e4tss100kb <- read.table(file = "data/factors-pmid38965376/annotated-rsid-p1e4/annotated_rsid_tss100kbp_p1e4.tsv",
                                          sep = "\t",
                                          header = TRUE) %>% 
  group_by(factor_id, gene_symbol) %>% 
  slice_min(order_by = pvalue, with_ties = FALSE) %>% 
  ungroup

factors_rsidGenes_P1e4locus <- read.table(file = "data/factors-pmid38965376/annotated-rsid-p1e4/annotated_rsid_locus0kbp_p1e4.tsv",
                                          sep = "\t",
                                          header = TRUE) %>% 
  group_by(factor_id, gene_symbol) %>% 
  slice_min(order_by = pvalue, with_ties = FALSE) %>% 
  ungroup

factors_rsidGenes_P1e4locus20kb <- read.table(file = "data/factors-pmid38965376/annotated-rsid-p1e4/annotated_rsid_locus20kbp_p1e4.tsv",
                                          sep = "\t",
                                          header = TRUE) %>% 
  group_by(factor_id, gene_symbol) %>% 
  slice_min(order_by = pvalue, with_ties = FALSE) %>% 
  ungroup

factors_rsidGenes_P1e4locus100kb <- read.table(file = "data/factors-pmid38965376/annotated-rsid-p1e4/annotated_rsid_locus100kbp_p1e4.tsv",
                                              sep = "\t",
                                              header = TRUE) %>% 
  group_by(factor_id, gene_symbol) %>% 
  slice_min(order_by = pvalue, with_ties = FALSE) %>% 
  ungroup

factors_rsidGenes_P1e4geneCenter100kb <- read.table(file = "data/factors-pmid38965376/annotated-rsid-p1e4/annotated_rsid_center100kbp_p1e4.tsv",
                                               sep = "\t",
                                               header = TRUE) %>% 
  group_by(factor_id, gene_symbol) %>% 
  slice_min(order_by = pvalue, with_ties = FALSE) %>% 
  ungroup
  

# ##############################################################################
# ---- analysis ----
# ##############################################################################
generate_and_permute_topn_factors_v2(
  gene_df         = lite_grSignatures,             # taka sama ramka genów jak wcześniej
  genes_vector    = hgnc_symbols_vector_v110,       # pełna lista symboli genów HGNC
  factor_data     = factors_rsidGenes_P1e4tss100kb, # Twoje dane z czynnikami
  n_randomization = 100,
  top_n           = 100,
  seed            = 1                            # opcjonalnie
) -> grSignature_factor_associations_top100P1e4tss100kb

generate_and_permute_topn_factors_v2(
  gene_df         = lite_grSignatures,             # taka sama ramka genów jak wcześniej
  genes_vector    = hgnc_symbols_vector_v110,       # pełna lista symboli genów HGNC
  factor_data     = factors_rsidGenes_P1e4locus, # Twoje dane z czynnikami
  n_randomization = 100,
  top_n           = 100,
  seed            = 1                            # opcjonalnie
) -> grSignature_factor_associations_top100P1e4locus

generate_and_permute_topn_factors_v2(
  gene_df         = lite_grSignatures,             # taka sama ramka genów jak wcześniej
  genes_vector    = hgnc_symbols_vector_v110,       # pełna lista symboli genów HGNC
  factor_data     = factors_rsidGenes_P1e4locus20kb, # Twoje dane z czynnikami
  n_randomization = 100,
  top_n           = 100,
  seed            = 1                            # opcjonalnie
) -> grSignature_factor_associations_top100P1e4locus20kb

generate_and_permute_topn_factors_v2(
  gene_df         = lite_grSignatures,             # taka sama ramka genów jak wcześniej
  genes_vector    = hgnc_symbols_vector_v110,       # pełna lista symboli genów HGNC
  factor_data     = factors_rsidGenes_P1e4locus100kb, # Twoje dane z czynnikami
  n_randomization = 100,
  top_n           = 100,
  seed            = 1                            # opcjonalnie
) -> grSignature_factor_associations_top100P1e4locus100kb


generate_and_permute_topn_factors_v2(
  gene_df         = lite_grSignatures,             # taka sama ramka genów jak wcześniej
  genes_vector    = hgnc_symbols_vector_v110,       # pełna lista symboli genów HGNC
  factor_data     = factors_rsidGenes_P1e4geneCenter100kb, # Twoje dane z czynnikami
  n_randomization = 100,
  top_n           = 100,
  seed            = 1                            # opcjonalnie
) -> grSignature_factor_associations_top100P1e4geneCenter100kb


# ##############################################################################
# ---- save to xlsx ----
# ##############################################################################
save_gr_signature_associations_to_xlsx(
  data_list   = grSignature_factor_associations_top100P1e4locus$original_association,
  output_file = "data/factors-pmid38965376/grSignaturesLite_top100P1e4random100locus.xlsx"
)

save_gr_signature_associations_to_xlsx(
  data_list   = grSignature_factor_associations_top100P1e4locus20kb$original_association,
  output_file = "data/factors-pmid38965376/grSignaturesLite_top100P1e4random100locus20kb.xlsx"
)

save_gr_signature_associations_to_xlsx(
  data_list   = grSignature_factor_associations_top100P1e4locus100kb$original_association,
  output_file = "data/factors-pmid38965376/grSignaturesLite_top100P1e4random100locus100kb.xlsx"
)

save_gr_signature_associations_to_xlsx(
  data_list   = grSignature_factor_associations_top100P1e4geneCenter100kb$original_association,
  output_file = "data/factors-pmid38965376/grSignaturesLite_top100P1e4random100geneCenter100kb.xlsx"
)


# ##############################################################################
# ---- summary results ----
# ##############################################################################

imap_dfr(
  list(
    tss100kb        = grSignature_factor_associations_top100P1e4tss100kb$original_association,
    locus           = grSignature_factor_associations_top100P1e4locus$original_association,
    locus20kb       = grSignature_factor_associations_top100P1e4locus20kb$original_association,
    locus100kb      = grSignature_factor_associations_top100P1e4locus100kb$original_association,
    geneCenter100kb = grSignature_factor_associations_top100P1e4geneCenter100kb$original_association
  ),
  function(signature_list, source_name) {
    imap_dfr(signature_list, function(sig_df, sig_name) {
      table(sig_df$fdr) %>%
        as.data.frame() %>%
        setNames(c("fdr", "freq")) %>%
        mutate(
          signature = sig_name,
          source = source_name
        )
    })
  }
) %>%
  complete(source, signature, fdr, fill = list(freq = 0)) %>%
  arrange(source, signature, fdr) %>% 
  pivot_wider(names_from = fdr, values_from = freq, values_fill = 0) -> all_fdr_summary_wide

all_fdr_summary_wide %>%
  pivot_longer(
    cols = starts_with("<") | starts_with("<="),
    names_to = "fdr",
    values_to = "freq"
  ) %>% 
  filter(fdr %in% c("< 0.05", "< 0.1", "< 0.2", "< 0.25", "< 0.3", "< 0.4")) %>% 
  ggplot(aes(x = source, y = freq, fill = signature)) +
  geom_bar(stat = "identity", position = "stack", width = 0.75) +
  facet_wrap(~ fdr, nrow = 2, strip.position = "top") +
  scale_fill_manual(
    values = c(
      "brain_down"         = "#FDBF6F",
      "brain_up"           = "#B2DF8A",
      "metasignature_down" = "#A6CEE3",
      "metasignature_up"   = "#CAB2D6"
    )
  ) +
  labs(
    title = "Associations by source (faceted by FDR threshold)",
    subtitle = "Each panel shows FDR threshold, bars are grouped by source and colored by signature",
    x = "Source of gene–variant mapping",
    y = "Number of associations",
    fill = "Signature"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(size = 18, hjust = 0.5, margin = margin(b = 15)),
    
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 18, face = "bold", margin = margin(r = 10)),
    
    axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 18),
    
    strip.text = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18),
    legend.position = "bottom"
  )
