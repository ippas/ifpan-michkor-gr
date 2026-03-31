
# ##############################################################################
# ---- chi2 analysis ----
# ##############################################################################
p1 <- process_overlap_heatmap_plot_only(factors_rsidGenes_P1e2locus20kb,
                                        type = "locus",
                                        window_kb = 20,
                                        plot_title_prefix = "A",
                                        pvalue_threshold = 0.000001,
                                        col_only_n_genes = TRUE,
                                        color_scale_range = c(0,5))

p2 <- process_overlap_heatmap_plot_only(factors_rsidGenes_P1e2locus30kb,
                                        type = "locus",
                                        window_kb = 30,
                                        plot_title_prefix = "B",
                                        pvalue_threshold = 0.000001,
                                        col_only_n_genes = TRUE,
                                        color_scale_range = c(0,5))


p3 <- process_overlap_heatmap_plot_only(factors_rsidGenes_P1e2locus50kb,
                                        type = "locus",
                                        window_kb = 50,
                                        plot_title_prefix = "C",
                                        pvalue_threshold = 0.000001,
                                        col_only_n_genes = TRUE,
                                        color_scale_range = c(0,5))

p4 <- process_overlap_heatmap_plot_only(factors_rsidGenes_P1e2locus75kb,
                                        type = "locus",
                                        window_kb = 75,
                                        plot_title_prefix = "D",
                                        pvalue_threshold = 0.000001,
                                        col_only_n_genes = TRUE,
                                        color_scale_range = c(0,5))

p5 <- process_overlap_heatmap_plot_only(factors_rsidGenes_P1e2locus100kb,
                                        type = "locus",
                                        window_kb = 100,
                                        plot_title_prefix = "E",
                                        pvalue_threshold = 0.000001,
                                        col_only_n_genes = TRUE,
                                        color_scale_range = c(0,5))


svg("data/factors-pmid38965376/figures/overlap-v1/combined_p1e6_min3_pthr_locus.svg", width = 28, height = 20)
(p1 / p2 / p3 / p4 / p5) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
dev.off()


p1 <- process_overlap_heatmap_plot_only(factors_rsidGenes_P1e2geneCenter20kb,
                                        type = "geneCenter",
                                        window_kb = 20,
                                        plot_title_prefix = "A",
                                        pvalue_threshold = 0.000001,
                                        col_only_n_genes = TRUE,
                                        color_scale_range = c(0,5))

p2 <- process_overlap_heatmap_plot_only(factors_rsidGenes_P1e2geneCenter30kb,
                                        type = "geneCenter",
                                        window_kb = 30,
                                        plot_title_prefix = "B",
                                        pvalue_threshold = 0.000001,
                                        col_only_n_genes = TRUE,
                                        color_scale_range = c(0,5))

p3 <- process_overlap_heatmap_plot_only(factors_rsidGenes_P1e2geneCenter50kb,
                                        type = "geneCenter",
                                        window_kb = 50,
                                        plot_title_prefix = "C",
                                        pvalue_threshold = 0.000001,
                                        col_only_n_genes = TRUE,
                                        color_scale_range = c(0,5))

p4 <- process_overlap_heatmap_plot_only(factors_rsidGenes_P1e2geneCenter75kb,
                                        type = "geneCenter",
                                        window_kb = 75,
                                        plot_title_prefix = "D",
                                        pvalue_threshold = 0.000001,
                                        col_only_n_genes = TRUE,
                                        color_scale_range = c(0,5))

p5 <- process_overlap_heatmap_plot_only(factors_rsidGenes_P1e2geneCenter100kb,
                                        type = "geneCenter",
                                        window_kb = 100,
                                        plot_title_prefix = "E",
                                        pvalue_threshold = 0.000001,
                                        col_only_n_genes = TRUE,
                                        color_scale_range = c(0,5))


svg("data/factors-pmid38965376/figures/overlap-v1/combined_p1e6_min3_pthr_center.svg", width = 28, height = 20)
(p1 / p2 / p3 / p4 / p5) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
dev.off()
