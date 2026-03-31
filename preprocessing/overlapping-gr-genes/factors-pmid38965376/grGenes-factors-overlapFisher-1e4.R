# --- Parametry ---
factors_df <- factors_rsidGenes_P1e2geneCenter100kb
type <- "geneCenter"
window_kb <- 100
plot_title_prefix <- "A"
pvalue_threshold <- 0.0001
drawing_overlap_threshold <- 3
fdr_output_threshold <- 0.2
col_only_n_genes <- TRUE
color_scale_range <- c(0, 4)

# --- Krok 1: przygotowanie listy genów ---
factor_gene_list <- factors_df %>%
  dplyr::filter(pvalue < pvalue_threshold) %>%
  dplyr::group_by(factor_name) %>%
  dplyr::summarise(genes = list(unique(gene_symbol)), .groups = "drop") %>%
  tibble::deframe()

gene_lists_all <- c(
  split(lite_grSignatures$hgnc_symbol, lite_grSignatures$signature_name),
  factor_gene_list
)

# --- Krok 2: test Fishera ---
fisher_results <- perform_fisher_tests(
  datasets = gene_lists_all,
  total_genes = hgnc_symbols_vector_v110
)

# --- Krok 3: przygotowanie filtrów ---
rows_to_filter <- unique(lite_grSignatures$signature_name)
cols_to_filter <- factors_df %>%
  dplyr::filter(pvalue < pvalue_threshold) %>%
  dplyr::pull(factor_name) %>%
  unique()

# --- Krok 4: przetwarzanie wyników ---
overlap <- processing_overlap_results_fisher(
  data = fisher_results,
  rows_to_filter = rows_to_filter,
  cols_to_filter = cols_to_filter,
  overlap_threshold = 0,
  fdr_threshold = 1,
  genes_list = gene_lists_all
)

draw_custom_heatmap_fisher_ggplot_v3(
  data_list = overlap,
  data_type = "original_data",
  p_thresholds = c(0.05, 0.01),
  color_rects = c("#4C8D05", "#66023C"),
  overlap_threshold = 4,
  title = "geneCenter +/- 100kb",
  col_only_n_genes = FALSE,
  apply_filling = F,
  palette = c("#2171b5", "#bdd7e7", "white", "#fcae91", "#cb181d"),
)


print_significant_rows(
  factors_df = factors_rsidGenes_P1e2geneCenter50kb,
  pvalue_threshold = 0.0001
)

