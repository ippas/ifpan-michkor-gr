large_gr_signatures <- list(
  gr_full_signatures_by_label = gene_list_log2ratio_top_50 %>%
    select(-top_50_genes) %>%
    unnest(data) %>%
    select(label_regulation, hgnc_symbol) %>%
    distinct() %>%
    group_by(label_regulation) %>%
    summarise(gene_list = list(unique(hgnc_symbol)), .groups = "drop") %>%
    deframe(),
  
  gr_full_signatures_by_tissue = gene_list_log2ratio_top_50 %>%
    select(-top_50_genes) %>%
    unnest(data) %>%
    ungroup() %>%
    select(simple_tissue, regulation, hgnc_symbol) %>%
    group_by(simple_tissue, regulation) %>%
    summarise(gene_list = list(unique(hgnc_symbol)), .groups = "drop") %>%
    mutate(name = paste(simple_tissue, regulation, sep = "_")) %>%
    select(name, gene_list) %>%
    deframe(),
  
  gr_top50_signatures_by_label = gene_list_log2ratio_top_50 %>%
    select(-data) %>%
    unnest(top_50_genes) %>%
    select(label_regulation, hgnc_symbol) %>%
    distinct() %>%
    group_by(label_regulation) %>%
    summarise(gene_list = list(unique(hgnc_symbol)), .groups = "drop") %>%
    deframe(),
  
  gr_top50_signatures_by_tissue = gene_list_log2ratio_top_50 %>%
    select(-data) %>%
    unnest(top_50_genes) %>%
    ungroup() %>%
    select(simple_tissue, regulation, hgnc_symbol) %>%
    group_by(simple_tissue, regulation) %>%
    summarise(gene_list = list(unique(hgnc_symbol)), .groups = "drop") %>%
    mutate(name = paste(simple_tissue, regulation, sep = "_")) %>%
    select(name, gene_list) %>%
    deframe()
)
