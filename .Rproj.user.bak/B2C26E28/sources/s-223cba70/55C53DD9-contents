


papers_data_preprocessing %>% 
  mutate(treatment = ifelse(treatment == "corticoterone", "corticosterone", treatment)) %>% 
  filter(!(treatment %in% c("NA", "PHA", "formoterol", "aldosterone", "vehicle-DMSO", "vitamin-d3", "mifepristone", "eplerenone"))) %>%
  filter((treatment %in% c("corticosterone", "cortisol", "prednisone", "prednisolone", "budesonide", "dexamethasone", "hydrocortisone"))) %>% 
  mutate(treatment_simple = ifelse(treatment %in% c("corticosterone", "cortisol"), "corticosterone_and_cortisol", treatment)) %>% 
  mutate(treatment_simple = ifelse(treatment %in% c("prednisone", "prednisolone"), "prednisone_and_prednisolone", treatment_simple)) %>%
  select(hgnc_symbol, treatment, treatment_simple) %>% 
  unique() -> venn_data

venn_data %>% 
  group_by(hgnc_symbol) %>% 
  nest() %>% 
  mutate(n_gc = map(data, ~nrow(.x))) %>% 
  unnest(n_gc) %>% 
  filter(n_gc == 5) %>% 
  .$hgnc_symbol -> common_gc_genes

venn_list <- venn_data %>%
  group_by(treatment_simple) %>%
  summarise(genes = list(unique(hgnc_symbol))) %>%
  mutate(gene_count = sapply(genes, length),
         treatment_label = paste0(treatment_simple, " (", gene_count, ")")) %>%
  select(treatment_label, genes) %>%
  deframe()

# Rysowanie wykresu Venna
venn.plot <- venn.diagram(
  x = venn_list, 
  filename = NULL, 
  fill = c("white", "white", "white", "white", "white"),
  alpha = 0.5,
  cat.col = c("black", "black", "black", "black", "black"),
  cat.cex = 1.2,
  main = "Venn Diagram - Glucocorticoids",
  main.cex = 1.5
)

# 40955 -> liczba recordów z bazy odnosząca się do tych 5 grup

venn_list %>% 
  lapply(., function(x) {x %>% unique() %>% length})

venn_list %>% unname() %>% unlist() %>% unique() %>% length()
# 12542 -> liczba genów z bazy odnosząca się do tych 5 grup


# Wyświetlanie wykresu
svg("results/figures/gene-database-summary/filtered-1-gene-database-summary//venn-treatment-gc.svg", width = 8, height = 8)
grid.draw(venn.plot)
dev.off()
