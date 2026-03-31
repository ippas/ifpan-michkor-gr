library(dplyr)
library(ggplot2)
library(tidyr)

# 1. Preliminary filtering and data preparation
filtered_data <- gene_list_log2ratio_top_50 %>% 
  select(-data) %>% 
  unnest(top_50_genes) %>% 
  ungroup %>% 
  mutate(log2ratio = as.numeric(log2ratio))


# 2. Create unique pairs of gene and regulation type
gene_regulation <- filtered_data %>%
  select(hgnc_symbol, regulation) %>%
  distinct()

# 3. Count the number of regulation types for each gene
gene_regulation_summary <- gene_regulation %>%
  count(hgnc_symbol, name = "n_regulation")

# 4. Extract genes that are regulated in both directions (bidirectional regulation)
bidirectionally_regulated_genes <- gene_regulation_summary %>%
  filter(n_regulation == 2) %>%
  pull(hgnc_symbol)

# 5. Extract genes regulated in only one direction:
#    - Only "down"
down_genes <- gene_regulation %>%
  group_by(hgnc_symbol) %>%
  filter(n() == 1, regulation == "down") %>%
  ungroup() %>%
  pull(hgnc_symbol)

#    - Only "up"
up_genes <- gene_regulation %>%
  group_by(hgnc_symbol) %>%
  filter(n() == 1, regulation == "up") %>%
  ungroup() %>%
  pull(hgnc_symbol)

# Display the number of genes in each category
cat("Number of bidirectionally regulated genes:", length(bidirectionally_regulated_genes), "\n")
cat("Number of genes regulated only 'up':", length(up_genes), "\n")
cat("Number of genes regulated only 'down':", length(down_genes), "\n")

# 6. Create a table showing counts of each regulation type for bidirectionally regulated genes
filtered_data %>%
  filter(hgnc_symbol %in% bidirectionally_regulated_genes) %>%
  count(regulation) %>%
  print()

# 7. Prepare data for visualization by adding a 'group' column
# Classify genes as:
#   - "up_only" / "down_only": if the gene is present in only one category,
#   - "up" / "down": if the gene is bidirectionally regulated
processed_gene_regulation <- filtered_data %>%
  mutate(group = case_when(
    hgnc_symbol %in% up_genes   ~ "up_only",
    hgnc_symbol %in% down_genes ~ "down_only",
    TRUE                        ~ regulation
  )) %>%
  mutate(group = factor(group, levels = c("down_only", "down", "up", "up_only")))

# 8. Visualize the data – create a boxplot of log2(FC) by gene regulation group
svg("/home/mateusz/projects/ifpan-michkor-gr/results/figures/harmonized-tissue-gene-lists/summary-gene-lists-top50/boxplot-up-down.svg", width = 6, height = 12)

ggplot(processed_gene_regulation, aes(x = group, y = log2ratio)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "Gene Regulation", y = "Log2(FC)")

dev.off()

# 1. Compute statistics per gene regulation group
# Compute statistics per gene regulation group
regulation_summary <- processed_gene_regulation %>%
  group_by(group) %>%
  dplyr::summarize(
    "Liczba wystąpień genów" = n(),                   # Total number of gene occurrences
    "Unikalne geny" = n_distinct(hgnc_symbol),        # Count of unique genes
    "Mediana log2ratio" = median(log2ratio, na.rm = TRUE) # Median log2ratio
  ) %>%
  pivot_longer(-group, names_to = "Metric", values_to = "Value") %>% # Convert to long format
  pivot_wider(names_from = group, values_from = Value) # Transpose to wide format

# Display results
print(regulation_summary)


