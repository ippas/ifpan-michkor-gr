# ##############################################################################
# ---- prepare data ----
# ##############################################################################
lite_grSignatures <- read.delim("results/gr-signatures/gr-signatures-multi-approach-24.04.2025.tsv",
           sep = "\t") %>% 
  filter(signature_name %in% c("universal_up", "universal_down", "brain_up", "brain_down")) %>% 
  mutate(signature_name = str_replace_all(signature_name, "universal", "metasignature"))

lite_grSignatures
hgnc_symbols_vector_v110

factors_rsidGenes_P1e2tss100kb <- read.table(file = "data/factors-pmid38965376/annotated-rsid-p1e2/annotated_rsid_tss100kbp_p1e2.tsv",
                                             sep = "\t",
                                             header = TRUE) %>% 
  group_by(factor_id, gene_symbol) %>% 
  slice_min(order_by = pvalue, with_ties = FALSE) %>% 
  ungroup

format(object.size(factors_rsidGenes_P1e2tss100kb), units = "auto")

factors_rsidGenes_P1e2locus <- read.table(file = "data/factors-pmid38965376/annotated-rsid-p1e2/annotated_rsid_locus0kbp_p1e2.tsv",
                                          sep = "\t",
                                          header = TRUE) %>% 
  group_by(factor_id, gene_symbol) %>% 
  slice_min(order_by = pvalue, with_ties = FALSE) %>% 
  ungroup

factors_rsidGenes_P1e2locus20kb <- read.table(file = "data/factors-pmid38965376/annotated-rsid-p1e2/annotated_rsid_locus20kbp_p1e2.tsv",
                                              sep = "\t",
                                              header = TRUE) %>% 
  group_by(factor_id, gene_symbol) %>% 
  slice_min(order_by = pvalue, with_ties = FALSE) %>% 
  ungroup


factors_rsidGenes_P1e2locus30kb <- read.table(file = "data/factors-pmid38965376/annotated-rsid-p1e2/annotated_rsid_locus30kbp_p1e2.tsv",
                                              sep = "\t",
                                              header = TRUE) %>% 
  group_by(factor_id, gene_symbol) %>% 
  slice_min(order_by = pvalue, with_ties = FALSE) %>% 
  ungroup

factors_rsidGenes_P1e2locus50kb <- read.table(file = "data/factors-pmid38965376/annotated-rsid-p1e2/annotated_rsid_locus50kbp_p1e2.tsv",
                                              sep = "\t",
                                              header = TRUE) %>% 
  group_by(factor_id, gene_symbol) %>% 
  slice_min(order_by = pvalue, with_ties = FALSE) %>% 
  ungroup


factors_rsidGenes_P1e2locus75kb <- read.table(file = "data/factors-pmid38965376/annotated-rsid-p1e2/annotated_rsid_locus75kbp_p1e2.tsv",
                                              sep = "\t",
                                              header = TRUE) %>% 
  group_by(factor_id, gene_symbol) %>% 
  slice_min(order_by = pvalue, with_ties = FALSE) %>% 
  ungroup


factors_rsidGenes_P1e2locus100kb <- read.table(file = "data/factors-pmid38965376/annotated-rsid-p1e2/annotated_rsid_locus100kbp_p1e2.tsv",
                                               sep = "\t",
                                               header = TRUE) %>% 
  group_by(factor_id, gene_symbol) %>% 
  slice_min(order_by = pvalue, with_ties = FALSE) %>% 
  ungroup



factors_rsidGenes_P1e2geneCenter20kb <- read.table(file = "data/factors-pmid38965376/annotated-rsid-p1e2/annotated_rsid_center20kbp_p1e2.tsv",
                                                   sep = "\t",
                                                   header = TRUE) %>% 
  group_by(factor_id, gene_symbol) %>% 
  slice_min(order_by = pvalue, with_ties = FALSE) %>% 
  ungroup


factors_rsidGenes_P1e2geneCenter30kb <- read.table(file = "data/factors-pmid38965376/annotated-rsid-p1e2/annotated_rsid_center30kbp_p1e2.tsv",
                                                    sep = "\t",
                                                    header = TRUE) %>% 
  group_by(factor_id, gene_symbol) %>% 
  slice_min(order_by = pvalue, with_ties = FALSE) %>% 
  ungroup


factors_rsidGenes_P1e2geneCenter50kb <- read.table(file = "data/factors-pmid38965376/annotated-rsid-p1e2/annotated_rsid_center50kbp_p1e2.tsv",
                                                   sep = "\t",
                                                   header = TRUE) %>% 
  group_by(factor_id, gene_symbol) %>% 
  slice_min(order_by = pvalue, with_ties = FALSE) %>% 
  ungroup


factors_rsidGenes_P1e2geneCenter75kb <- read.table(file = "data/factors-pmid38965376/annotated-rsid-p1e2/annotated_rsid_center75kbp_p1e2.tsv",
                                                   sep = "\t",
                                                   header = TRUE) %>% 
  group_by(factor_id, gene_symbol) %>% 
  slice_min(order_by = pvalue, with_ties = FALSE) %>% 
  ungroup

factors_rsidGenes_P1e2geneCenter100kb <- read.table(file = "data/factors-pmid38965376/annotated-rsid-p1e2/annotated_rsid_center100kbp_p1e2.tsv",
                                                    sep = "\t",
                                                    header = TRUE) %>% 
  group_by(factor_id, gene_symbol) %>% 
  slice_min(order_by = pvalue, with_ties = FALSE) %>% 
  ungroup


factors_rsidGenes_P1e2tss100kb %>% filter(pvalue < 0.0001) %>% dim
factors_rsidGenes_P1e2locus %>% filter(pvalue < 0.0001) %>% dim
factors_rsidGenes_P1e2locus20kb %>% filter(pvalue < 0.0001) %>% dim
factors_rsidGenes_P1e2locus100kb %>% filter(pvalue < 0.0001) %>% dim
factors_rsidGenes_P1e2geneCenter100kb%>% filter(pvalue < 0.0001) %>% dim

# ##############################################################################
# ---- analysis ----
# ##############################################################################
generate_and_permute_topn_factors_v2(
  gene_df         = lite_grSignatures,             # taka sama ramka genów jak wcześniej
  genes_vector    = hgnc_symbols_vector_v110,       # pełna lista symboli genów HGNC
  factor_data     = factors_rsidGenes_P1e2tss100kb, # Twoje dane z czynnikami
  n_randomization = 1000,
  top_n           = 400,
  seed            = 1                            # opcjonalnie
) -> grSignature_factor_associations_top400P1e2tss100kb

grSignature_factor_associations_top400P1e2tss100kb$original_association %>% 
  lapply(., function(x) {
    x %>% filter(pvalue < 0.0001) %>% dim
  })

generate_and_permute_topn_factors_v2(
  gene_df         = lite_grSignatures,             # taka sama ramka genów jak wcześniej
  genes_vector    = hgnc_symbols_vector_v110,       # pełna lista symboli genów HGNC
  factor_data     = factors_rsidGenes_P1e2locus, # Twoje dane z czynnikami
  n_randomization = 1000,
  top_n           = 400,
  seed            = 1                            # opcjonalnie
) -> grSignature_factor_associations_top400P1e2locus

grSignature_factor_associations_top400P1e2locus$original_association %>% 
  lapply(., function(x) {
    x %>% filter(pvalue < 0.0001) %>% dim
  })

generate_and_permute_topn_factors_v2(
  gene_df         = lite_grSignatures,             # taka sama ramka genów jak wcześniej
  genes_vector    = hgnc_symbols_vector_v110,       # pełna lista symboli genów HGNC
  factor_data     = factors_rsidGenes_P1e2locus20kb, # Twoje dane z czynnikami
  n_randomization = 1000,
  top_n           = 400,
  seed            = 1                            # opcjonalnie
) -> grSignature_factor_associations_top400P1e2locus20kb

grSignature_factor_associations_top400P1e2locus20kb$original_association %>% 
  lapply(., function(x) {
    x %>% filter(pvalue < 0.0001) %>% dim
  })

generate_and_permute_topn_factors_v2(
  gene_df         = lite_grSignatures,             # taka sama ramka genów jak wcześniej
  genes_vector    = hgnc_symbols_vector_v110,       # pełna lista symboli genów HGNC
  factor_data     = factors_rsidGenes_P1e2locus100kb, # Twoje dane z czynnikami
  n_randomization = 1000,
  top_n           = 400,
  seed            = 1                            # opcjonalnie
) -> grSignature_factor_associations_top400P1e2locus100kb

grSignature_factor_associations_top400P1e2locus100kb$original_association %>% 
  lapply(., function(x) {
    x %>% filter(pvalue < 0.0001) %>% dim
  })

generate_and_permute_topn_factors_v2(
  gene_df         = lite_grSignatures,             # taka sama ramka genów jak wcześniej
  genes_vector    = hgnc_symbols_vector_v110,       # pełna lista symboli genów HGNC
  factor_data     = factors_rsidGenes_P1e2geneCenter100kb, # Twoje dane z czynnikami
  n_randomization = 1000,
  top_n           = 400,
  seed            = 1                            # opcjonalnie
) -> grSignature_factor_associations_top400P1e2geneCenter100kb

grSignature_factor_associations_top400P1e2geneCenter100kb$original_association %>% 
  lapply(., function(x) {
    x %>% filter(pvalue < 0.0001)
  })


z# rm(
#   grSignature_factor_associations_top100P1e2tss100kb,
#   grSignature_factor_associations_top100P1e2locus,
#   grSignature_factor_associations_top100P1e2locus20kb,
#   grSignature_factor_associations_top100P1e2locus100kb,
#   grSignature_factor_associations_top100P1e2geneCenter100kb
# )
# ##############################################################################
# ---- save to xlsx ----
# ##############################################################################
save_gr_signature_associations_to_xlsx(
  data_list   = grSignature_factor_associations_top400P1e2tss100kb$original_association %>% 
    lapply(., function(x){ x %>% filter(pvalue < 0.0001)}),
  output_file = "data/factors-pmid38965376/annotated-rsid-p1e2/grSignaturesLite_P1e4random1000tss100kb.xlsx"
)

save_gr_signature_associations_to_xlsx(
  data_list   = grSignature_factor_associations_top400P1e2locus$original_association %>% 
    lapply(., function(x){ x %>% filter(pvalue < 0.0001)}),
  output_file = "data/factors-pmid38965376/annotated-rsid-p1e2/grSignaturesLite_P1e4random1000locus.xlsx"
)

save_gr_signature_associations_to_xlsx(
  data_list   = grSignature_factor_associations_top400P1e2locus20kb$original_association %>% 
    lapply(., function(x){ x %>% filter(pvalue < 0.0001)}),
  output_file = "data/factors-pmid38965376/annotated-rsid-p1e2/grSignaturesLite_P1e4random1000locus20kb.xlsx"
)

save_gr_signature_associations_to_xlsx(
  data_list   = grSignature_factor_associations_top400P1e2locus100kb$original_association %>% 
    lapply(., function(x){ x %>% filter(pvalue < 0.0001)}),
  output_file = "data/factors-pmid38965376/annotated-rsid-p1e2/grSignaturesLite_P1e4random1000locus100kb.xlsx"
)

save_gr_signature_associations_to_xlsx(
  data_list   = grSignature_factor_associations_top400P1e2geneCenter100kb$original_association %>% 
    lapply(., function(x){ x %>% filter(pvalue < 0.0001)}),
  output_file = "data/factors-pmid38965376/annotated-rsid-p1e2/grSignaturesLite_P1e4random1000geneCenter100kb.xlsx"
)


# ##############################################################################
# ---- summary results ----
# ##############################################################################

imap_dfr(
  list(
    tss100kb        = grSignature_factor_associations_top400P1e2tss100kb$original_association %>% lapply(., function(x){ x %>% filter(pvalue < 0.0001)}),
    locus           = grSignature_factor_associations_top400P1e2locus$original_association %>% lapply(., function(x){ x %>% filter(pvalue < 0.0001)}),
    locus20kb       = grSignature_factor_associations_top400P1e2locus20kb$original_association %>% lapply(., function(x){ x %>% filter(pvalue < 0.0001)}),
    locus100kb      = grSignature_factor_associations_top400P1e2locus100kb$original_association %>% lapply(., function(x){ x %>% filter(pvalue < 0.0001)}),
    geneCenter100kb = grSignature_factor_associations_top400P1e2geneCenter100kb$original_association %>% lapply(., function(x){ x %>% filter(pvalue < 0.0001)})
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


all_fdr_summary_wide %>%
  pivot_longer(
    cols = starts_with("<") | starts_with("<="),
    names_to = "fdr_raw",
    values_to = "freq"
  ) %>%
  mutate(fdr_numeric = readr::parse_number(fdr_raw)) %>%
  arrange(source, signature, fdr_numeric) %>%
  group_by(source, signature) %>%
  mutate(freq_cumsum = cumsum(freq)) %>%
  ungroup() %>%
  mutate(fdr_bin = case_when(
    fdr_numeric < 0.05 ~ "FDR < 0.05",
    fdr_numeric < 0.1  ~ "0.05 ≤ FDR < 0.1",
    fdr_numeric < 0.2  ~ "0.1 ≤ FDR < 0.2",
    fdr_numeric < 0.3  ~ "0.2 ≤ FDR < 0.3",
    fdr_numeric < 0.4  ~ "0.3 ≤ FDR < 0.4",
    TRUE               ~ "FDR ≥ 0.4"
  )) %>%
  group_by(source, signature, fdr_bin) %>%
  summarise(freq = max(freq_cumsum), .groups = "drop") %>%
  ggplot(aes(x = source, y = freq, fill = signature)) +
  geom_bar(stat = "identity", position = "stack", width = 0.75) +
  facet_wrap(~ fdr_bin, nrow = 2, strip.position = "top") +
  scale_fill_manual(
    values = c(
      "brain_down"         = "#FDBF6F",
      "brain_up"           = "#B2DF8A",
      "metasignature_down" = "#A6CEE3",
      "metasignature_up"   = "#CAB2D6"
    )
  ) +
  labs(
    title = "Associations by source (faceted by FDR bins)",
    subtitle = "Each panel shows FDR interval, bars are grouped by source and colored by signature",
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
