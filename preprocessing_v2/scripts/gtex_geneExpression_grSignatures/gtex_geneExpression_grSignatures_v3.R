
hgnc_symbols_vector_v110_gtex <- multi_downloadMedianSubtissuesGTEx_tissueSummaryExpression(gene_symbol = hgnc_symbols_vector_v110)


write_tsv(
  hgnc_symbols_vector_v110_gtex,
  "/home/mateusz/projects/ifpan-michkor-gr/data/gtex/proteinCodingGenes_biomartv110_gtexExpression_25.02.2026.tsv"
)

grSystemicTissues_all %>% 
  filter(grSignature == "lungDown") %>% 
  filter(p_value < 0.05) %>% 
  filter(observed_overlap >= 3) %>% 
  .$overlap_genes %>% 
  strsplit(",") %>% 
  unlist() %>% 
  unique -> grLungDown_genesSignif

grLungDowGenesSignif_gtex <- multi_downloadMedianSubtissuesGTEx_tissueSummaryExpression(gene_symbol = grLungDown_genesSignif)

calculate_tau <- function(x, min_expr = 1) {
  
  # usuń geny praktycznie niewyrażane
  if (max(x, na.rm = TRUE) < min_expr) return(NA_real_)
  
  x_norm <- x / max(x, na.rm = TRUE)
  
  tau <- sum(1 - x_norm, na.rm = TRUE) / (length(x_norm) - 1)
  
  return(tau)
}

grLungDowGenesSignif_gtex %>% 
  select(geneSymbol, tissue, median_max) %>% 
  group_by(geneSymbol) %>% 
  nest() %>% 
  mutate(expression = map(data, ~ .x$median_max)) %>% 
  mutate(tau_thr1 = map(data, ~ calculate_tau(.x$median_max, min_expr = 50))) %>% .$tau_thr1 %>% unlist %>% summary


grSystemicTissues_all %>% 
  filter(grSignature == "neuralUp") %>% 
  filter(p_value < 0.05) %>% 
  filter(observed_overlap >= 3) %>% 
  .$overlap_genes %>% 
  strsplit(",") %>% 
  unlist() %>% 
  unique -> grNeuralUp_genesSignif

grNeuralUpGenesSignif_gtex <- multi_downloadMedianSubtissuesGTEx_tissueSummaryExpression(gene_symbol = grNeuralUp_genesSignif)

grNeuralUpGenesSignif_gtex %>% 
  select(geneSymbol, tissue, median_max) %>% 
  group_by(geneSymbol) %>% 
  nest() %>% 
  mutate(expression = map(data, ~ .x$median_max)) %>% 
  mutate(tau_thr1 = map(data, ~ calculate_tau(.x$median_max))) %>% .$tau_thr1 %>% unlist %>% summary
