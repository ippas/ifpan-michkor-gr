# ##############################################################################
# ---- prepare data ----
# ##############################################################################
factors_rsidGenes_P1e4tss100kb <- read.table(file = "data/factors-pmid38965376/factorID-pmid38965376-rsidP1e4-proteinCodingv110-tss100kb.tsv", 
                                             sep = "\t", 
                                             header = TRUE) 


factors_rsidGenes_P1e4tss100kb %>% 
  group_by(factor_id, gene_symbol) %>% 
  slice_min(order_by = pvalue, with_ties = FALSE) %>% 
  ungroup-> factors_rsidGenes_P1e4tss100kb_minP



lite_grSignatures
hgnc_symbols_vector_v110


# ##############################################################################
# ---- functions ----
# ##############################################################################
source("preprocessing/overlapping-gr-genes/factors-pmid38965376/factorsPMID38965376-functions.R")

# ##############################################################################
# ---- analysis ----
# ##############################################################################
generate_and_permute_topn_factors_v2(
  gene_df         = lite_grSignatures,             # taka sama ramka genów jak wcześniej
  genes_vector    = hgnc_symbols_vector_v110,       # pełna lista symboli genów HGNC
  factor_data     = factors_rsidGenes_P1e4tss100kb_minP, # Twoje dane z czynnikami
  n_randomization = 100,
  top_n           = 200,
  seed            = 1                            # opcjonalnie
) -> grSignature_factor_associations_FDRMonteCarlo


grSignature_factor_associations_FDRMonteCarlo$original_association[
  c("metasignature_up", "metasignature_down", "brain_up", "brain_down")
] <- grSignature_factor_associations_FDRMonteCarlo$original_association[
  c("metasignature_up", "metasignature_down", "brain_up", "brain_down")
] %>%
  map(~ .x %>%
        as.data.frame() %>%
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
        ))



grSignature_factor_associations_FDRMonteCarlo$original_association$brain_down %>% 
  filter(fdr %in% c("< 0.05", "< 0.1","< 0.2")) %>% .$gene_symbol %>% unique %>% length()
# ##############################################################################
# ---- save to xlsx ----
# ##############################################################################

save_gr_signature_associations_to_xlsx(
  data_list = grSignature_factor_associations_FDRMonteCarlo$original_association,
  output_file = "data/factors-pmid38965376/grSignaturesLite_factorsAssocitionP1e4_withRandomFDRtop200random100-07.07.2025.xlsx"
)
