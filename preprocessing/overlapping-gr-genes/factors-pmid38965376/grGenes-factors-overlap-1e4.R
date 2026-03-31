# ##############################################################################
# ---- prepare data ----
# ##############################################################################

lite_grSignatures
hgnc_symbols_vector_v110

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



papers_data_preprocessing %>% filter(simple_tissue == "brain") %>% 
  filter(regulation == "up") %>%
  .$hgnc_symbol %>% table() %>% 
  as.data.frame() %>%  
  set_colnames(c("hgnc_symbol", "freq")) %>% 
  filter(freq > 2) %>% 
  mutate(signature_name = "brain_up_minFreq3") %>% 
  mutate(signature_derivation = "tissue_freq") %>% 
  select(-freq) -> brain_up_minFreq3

papers_data_preprocessing %>% filter(simple_tissue == "brain") %>% 
  filter(regulation == "down") %>%
  .$hgnc_symbol %>% table() %>% 
  as.data.frame() %>%  
  set_colnames(c("hgnc_symbol", "freq")) %>% 
  filter(freq > 2) %>% 
  mutate(signature_name = "brain_down_minFreq3") %>% 
  mutate(signature_derivation = "tissue_freq") %>% 
  select(-freq) -> brain_down_minFreq3

lite_grSignatures %>% 
  filter(signature_derivation != "tissue_freq") -> lite_grSignatures

rbind(lite_grSignatures, brain_up_minFreq3, brain_down_minFreq3) -> lite_grSignatures


# ##############################################################################
# ---- chi2 analysis ----
# ##############################################################################
p1 <- process_overlap_heatmap_plot_only(factors_rsidGenes_P1e2locus20kb,
                                        type = "locus",
                                        window_kb = 20,
                                        plot_title_prefix = "A",
                                        pvalue_threshold = 0.0001,
                                        color_scale_range = c(0, 8),
                                        col_only_n_genes = TRUE)

p2 <- process_overlap_heatmap_plot_only(factors_rsidGenes_P1e2locus30kb,
                                        type = "locus",
                                        window_kb = 30,
                                        plot_title_prefix = "B",
                                        pvalue_threshold = 0.0001,
                                        col_only_n_genes = TRUE)


p3 <- process_overlap_heatmap_plot_only(factors_rsidGenes_P1e2locus50kb,
                                        type = "locus",
                                        window_kb = 50,
                                        plot_title_prefix = "C",
                                        pvalue_threshold = 0.0001,
                                        col_only_n_genes = TRUE)

p4 <- process_overlap_heatmap_plot_only(factors_rsidGenes_P1e2locus75kb,
                                        type = "locus",
                                        window_kb = 75,
                                        plot_title_prefix = "D",
                                        pvalue_threshold = 0.0001,
                                        col_only_n_genes = TRUE)

p5 <- process_overlap_heatmap_plot_only(factors_rsidGenes_P1e2locus100kb,
                                        type = "locus",
                                        window_kb = 100,
                                        plot_title_prefix = "E",
                                        pvalue_threshold = 0.0001,
                                        col_only_n_genes = TRUE)


svg("data/factors-pmid38965376/figures/overlap-v1/combined_p1e4_min3_pthr_locus.svg", width = 28, height = 20)
(p1 / p2 / p3 / p4 / p5) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
dev.off()


p1 <- process_overlap_heatmap_plot_only(factors_rsidGenes_P1e2geneCenter20kb,
                                        type = "geneCenter",
                                        window_kb = 20,
                                        plot_title_prefix = "A",
                                        pvalue_threshold = 0.0001,
                                        col_only_n_genes = TRUE)

p2 <- process_overlap_heatmap_plot_only(factors_rsidGenes_P1e2geneCenter30kb,
                                        type = "geneCenter",
                                        window_kb = 30,
                                        plot_title_prefix = "B",
                                        pvalue_threshold = 0.0001,
                                        col_only_n_genes = TRUE)

p3 <- process_overlap_heatmap_plot_only(factors_rsidGenes_P1e2geneCenter50kb,
                                        type = "geneCenter",
                                        window_kb = 50,
                                        plot_title_prefix = "C",
                                        pvalue_threshold = 0.0001,
                                        col_only_n_genes = TRUE)

p4 <- process_overlap_heatmap_plot_only(factors_rsidGenes_P1e2geneCenter75kb,
                                        type = "geneCenter",
                                        window_kb = 75,
                                        plot_title_prefix = "D",
                                        pvalue_threshold = 0.0001,
                                        col_only_n_genes = TRUE)

process_overlap_heatmap_plot_only(factors_rsidGenes_P1e2geneCenter75kb,
                                  type = "geneCenter",
                                  window_kb = 75,
                                  plot_title_prefix = "D",
                                  pvalue_threshold = 0.0001,
                                  col_only_n_genes = TRUE, return_df = T) %>% 
  .$df %>% 
  filter(Var1 %in% c("brain_up", "brain_down")) %>% 
  filter(Var2 %in% c("alcohol_use_and_misuse", "anxiety_and_nervousness",
                     "clinical_anxiety_and_depression", "cognition_and_processing_speed",
                     "depressive_symptomatology", "trauma", "social_and_economic_stability")) %>% 
  filter(gene_overlap_count > 2)

p5 <- process_overlap_heatmap_plot_only(factors_rsidGenes_P1e2geneCenter100kb,
                                        type = "geneCenter",
                                        window_kb = 100,
                                        plot_title_prefix = "E",
                                        pvalue_threshold = 0.001,
                                        col_only_n_genes = TRUE)



process_overlap_heatmap_plot_only(factors_rsidGenes_P1e2geneCenter100kb,
                                  type = "geneCenter",
                                  window_kb = 100,
                                  plot_title_prefix = "E",
                                  pvalue_threshold = 0.0001,
                                  col_only_n_genes = F,
                                  return_df = T) %>% 
  .$df %>% 
  filter(Var1 %in% c("brain_up", "brain_down")) %>% 
  filter(Var2 %in% c("alcohol_use_and_misuse", "anxiety_and_nervousness",
                     "clinical_anxiety_and_depression", "cognition_and_processing_speed",
                     "depressive_symptomatology", "trauma", "social_and_economic_stability")) %>% 
  filter(gene_overlap_count > 2)



svg("data/factors-pmid38965376/figures/overlap-v1/combined_p1e4_min3_pthr_center.svg", width = 28, height = 20)
(p1 / p2 / p3 / p4 / p5) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
dev.off()


