# ##############################################################################
# ---- prepare data ----
# ##############################################################################
read.delim(gzfile("data/PGC/mergedPGC_annotation_geneCenter50kb_p1e4.tsv"),
           header = TRUE, sep = "\t") %>% 
  filter(source_file != "") -> pgc_annotation_geneCenter50kb_p1e4

pgc_annotation_geneCenter50kb_p1e4 %>% .$gene_symbol %>% unique() %>% length()

pgc_annotation_geneCenter50kb_p1e4 %>% 
  filter(pvalue < 0.00000001) %>% 
  # head(1000) %>% 
  select( source_file, gene_symbol) %>% 
  unique %>% 
  group_by(source_file) %>% 
  nest() %>% 
  mutate(n_genes = map(data, ~ .x %>% nrow)) %>% 
  unnest() %>% 
  .$n_genes %>% sd



AllBrain2BrainSignatures2Global %>% length()

pgcGeneCenter50kb_AllBrain2Brain2Global_rsIDp1e4 <- process_overlap_heatmap_plot_only(pgc_annotation_geneCenter50kb_p1e4,
                                                                          type = "geneCenter",
                                                                          window_kb = 50,
                                                                          reference_gene_lists = AllBrain2BrainSignatures2Global,
                                                                          plot_title_prefix = "A",
                                                                          pvalue_threshold = 0.0001,
                                                                          color_scale_range = c(0, 8),
                                                                          phenotype_col = "source_file",
                                                                          verbose = T,
                                                                          col_only_n_genes = F)

pgcGeneCenter50kb_AllBrain2Brain2Global_rsIDp1e4$overlap$original_data$df %>% 
  filter(p_value < 0.001) %>% 
  filter(gene_overlap_count > 2) %>% 
  filter(!(Var1 %in% c("brain_down", "brain_up", "metasignature_down", "metasignature_up"))) %>% 
  .$overlap_genes %>% 
  strsplit(",") %>% 
  unlist %>% 
  table %>% 
  as.data.frame() %>% 
  arrange(desc(Freq)) %>% 
  set_colnames(c("gene_symbol", "freq")) %>% 
  filter(gene_symbol %in% {read.delim(file = "results/gr-signatures/gr-signatures-multi-approach-24.04.2025.tsv") %>% 
           filter(signature_derivation == "marpiech_cluster") %>% 
           select(-signature_derivation) %>% 
           filter(signature_name == "cluster_O") %>% .$hgnc_symbol})

pgcGeneCenter50kb_AllBrain2Brain2Global_rsIDp1e5 <- process_overlap_heatmap_plot_only(pgc_annotation_geneCenter50kb_p1e4,
                                                                                      type = "geneCenter",
                                                                                      window_kb = 50,
                                                                                      reference_gene_lists = AllBrain2BrainSignatures2Global,
                                                                                      plot_title_prefix = "A",
                                                                                      pvalue_threshold = 0.00001,
                                                                                      color_scale_range = c(0, 8),
                                                                                      phenotype_col = "source_file",
                                                                                      verbose = T,
                                                                                      col_only_n_genes = F)

pgcGeneCenter50kb_AllBrain2Brain2Global_rsIDp1e6 <- process_overlap_heatmap_plot_only(pgc_annotation_geneCenter50kb_p1e4,
                                                                                      type = "geneCenter",
                                                                                      window_kb = 50,
                                                                                      reference_gene_lists = AllBrain2BrainSignatures2Global,
                                                                                      plot_title_prefix = "A",
                                                                                      pvalue_threshold = 0.000001,
                                                                                      color_scale_range = c(0, 8),
                                                                                      phenotype_col = "source_file",
                                                                                      verbose = T,
                                                                                      col_only_n_genes = F)

pgcGeneCenter50kb_AllBrain2Brain2Global_rsIDp1e7 <- process_overlap_heatmap_plot_only(pgc_annotation_geneCenter50kb_p1e4,
                                                                                      type = "geneCenter",
                                                                                      window_kb = 50,
                                                                                      reference_gene_lists = AllBrain2BrainSignatures2Global,
                                                                                      plot_title_prefix = "A",
                                                                                      pvalue_threshold = 0.0000001,
                                                                                      color_scale_range = c(0, 8),
                                                                                      phenotype_col = "source_file",
                                                                                      verbose = T,
                                                                                      col_only_n_genes = F)


pgcGeneCenter50kb_AllBrain2Brain2Global_rsIDp1e8 <- process_overlap_heatmap_plot_only(pgc_annotation_geneCenter50kb_p1e4,
                                                                                      type = "geneCenter",
                                                                                      window_kb = 50,
                                                                                      reference_gene_lists = AllBrain2BrainSignatures2Global,
                                                                                      plot_title_prefix = "A",
                                                                                      pvalue_threshold = 0.00000001,
                                                                                      color_scale_range = c(0, 8),
                                                                                      phenotype_col = "source_file",
                                                                                      verbose = T,
                                                                                      col_only_n_genes = F)


pgcGeneCenter50kb_AllBrain2Brain2Global_rsIDp1e8$overlap$original_data$df %>% 
  filter(p_value < 0.01, gene_overlap_count  > 2)


AllBrain2BrainSignatures2GlobalMarpiechClusters[[c("cluster_P")]]

pgcGeneCenter50kb_AllBrain2Brain2GlobalMarpiechClusters_rsIDp1e4 <- process_overlap_heatmap_plot_only(pgc_annotation_geneCenter50kb_p1e4,
                                                                                      type = "geneCenter",
                                                                                      window_kb = 50,
                                                                                      reference_gene_lists = AllBrain2BrainSignatures2GlobalMarpiechClusters,
                                                                                      plot_title_prefix = "A",
                                                                                      pvalue_threshold = 0.0001,
                                                                                      color_scale_range = c(0, 8),
                                                                                      phenotype_col = "source_file",
                                                                                      verbose = T,
                                                                                      col_only_n_genes = F)


pgcGeneCenter50kb_AllBrain2Brain2GlobalMarpiechClusters_rsIDp1e4$overlap$original_data$df %>% 
  filter(Var1 %in% c("cluster_A", "cluster_B", "cluster_C", "cluster_D", "cluster_E", "cluster_F", 
                     "cluster_G", "cluster_H", "cluster_I", "cluster_J", "cluster_K", "cluster_L", 
                     "cluster_M", "cluster_N", "cluster_O", "cluster_P")) %>% 
  # filter(Var1 %in% c("cluster_O", "cluster_P", "cluster_K")) %>% 
  filter(p_value < 0.01, gene_overlap_count > 2)
