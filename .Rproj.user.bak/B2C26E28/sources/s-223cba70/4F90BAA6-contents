# 1. Filter and prepare the d
venn_data <- filt_1_12_gr_database %>%
  # Keep only records with regulation "up" or "down"
  filter(regulation %in% c("up", "down")) %>%
  # Exclude unwanted treatments
  filter(!(treatment %in% c("NA", "PHA", "formoterol", "aldosterone", 
                            "vehicle-DMSO", "vitamin-d3", "mifepristone", "eplerenone"))) %>%
  # Keep only records with treatments corresponding to glucocorticoids
  filter(treatment %in% c("corticosterone", "cortisol", "prednisone", 
                          "prednisolone", "budesonide", "dexamethasone", "hydrocortisone")) %>%
  # Simplify treatment names using case_when
  mutate(
    treatment_simple = case_when(
      treatment %in% c("corticosterone", "cortisol", "hydrocortisone") ~ "corticosterone_and_cortisol",
      treatment %in% c("prednisone", "prednisolone") ~ "prednisone_and_prednisolone",
      TRUE ~ treatment
    )
  ) %>%
  # Select relevant columns and remove duplicate rows
  select(hgnc_symbol, treatment, treatment_simple) %>%
  distinct()

# 2. Identify genes present in all 5 glucocorticoid groups
common_gc_genes <- venn_data %>%
  group_by(hgnc_symbol) %>%
  nest() %>%
  mutate(n_gc = map_int(data, nrow)) %>%
  ungroup() %>%
  filter(n_gc == 5) %>%
  pull(hgnc_symbol)

# 3. Create a list for the Venn diagram with gene counts in the names
venn_list <- venn_data %>%
  group_by(treatment_simple) %>%
  summarise(genes = list(unique(hgnc_symbol)), .groups = "drop") %>%
  # Calculate the number of genes in each treatment group
  mutate(count = map_int(genes, length)) %>%
  # Create new labels that include the gene count
  mutate(treatment_label = paste0(treatment_simple, " (", count, ")")) %>%
  select(treatment_label, genes) %>%
  deframe()

# Print the number of unique genes in each treatment group
print(lapply(venn_list, function(x) length(unique(x))))

# Calculate the total number of unique genes across the 5 groups
total_genes <- venn_list %>%
  unname() %>%
  unlist() %>%
  unique() %>%
  length()
print(total_genes)

# 4. Draw the Venn diagram
venn_plot <- venn.diagram(
  x = venn_list,
  filename = NULL,
  fill = rep("white", length(venn_list)),
  alpha = 0.5,
  cat.col = rep("black", length(venn_list)),
  cat.cex = 1.2,
  main = "Venn Diagram - Glucocorticoids",
  main.cex = 1.5
)

# Create a new page and draw the Venn diagram
grid.newpage()
grid.draw(venn_plot)

dev.off()
# Optionally: Save the plot as an SVG file
svg("results/figures/gene-database-summary/filtered-1-12-gene-database-summary//venn-treatment-gc.svg", width = 8, height = 8)
grid.draw(venn_plot)
dev.off()


venn_list %>%
  unname() %>% 
  unlist %>% 
  table() %>% 
  as.data.frame() %>% 
  set_colnames(c("hgnc_symbol", "n_gc")) %>% 
  filter(n_gc == 4)


stack(venn_list) %>% 
  set_colnames(c("hgnc_symbol", "gc")) %>% 
  group_by(hgnc_symbol) %>% 
  nest %>% 
  mutate(n_gc = map(data, ~nrow(.x))) %>% 
  unnest(n_gc) %>% 
  filter(n_gc == 3) %>% 
  unnest %>% 
  filter(gc != "dexamethasone (663)") %>% 
  select(-n_gc) %>% 
  group_by(hgnc_symbol) %>% 
  nest %>% 
  mutate(n_gc = map(data, ~nrow(.x))) %>% 
  unnest(n_gc) %>% filter(n_gc ==3)
  


venn_list %>%
  unname() %>% 
  unlist %>% 
  table() %>% 
  as.data.frame() %>% 
  set_colnames(c("hgnc_symbol", "n_gc")) %>% 
  filter(n_gc == 1) %>% 
  .$hgnc_symbol
