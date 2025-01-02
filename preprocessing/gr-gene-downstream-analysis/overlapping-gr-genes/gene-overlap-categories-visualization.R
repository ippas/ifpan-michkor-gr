psychiatric_association_category <- c("Mental distress", "Mental health", 'Smoking', "Psychosocial factors", "Anxiety", "Depression", "Alcohol")

overlap_resutls_random_dependent_fdr0.01 %>% 
  # Convert the list to a tibble with each list element as a row
  enframe(name = "category", value = "data") %>%  
  mutate(signif_overlap = map(data, ~nrow(.x$significant_uniq_data$df))) %>% 
  unnest(signif_overlap) %>% 
  mutate(permutation_FDR = map2(data, signif_overlap, ~sum(.x$permutation_results >= .y)/1000)) %>% 
  unnest(permutation_FDR) %>% 
  mutate(signif_phenotypes = map(data, ~ .x$significant_uniq_data$df$Var1 %>%  unique() %>% length())) %>% 
  mutate(signif_papers = map(data, ~ .x$significant_uniq_data$df$Var2 %>%  unique() %>% length())) %>% 
  mutate(number_phenotypes = map(data, ~ .x$original_data$cols %>% length())) %>% 
  unnest(c(signif_phenotypes, signif_papers, number_phenotypes)) %>% 
  mutate(ratio_signif_phenotypes = signif_phenotypes/number_phenotypes) %>% 
  mutate(FDR_threshold = ifelse(permutation_FDR > 0.05, ">0.05", 
                                ifelse(permutation_FDR < 0.05 & permutation_FDR > 0.001, "<0.05 and >0.001", "<0.001"))) %>% 
  mutate(colors_psychiatric = ifelse(category %in% psychiatric_association_category, "#6495ED", "black"))  -> overlap_results_preprocessing


# Plot 1: Visualize the number of significant overlaps for each category
plot1 <-
  gene_overlap_barplot_summary(
    overlap_results_preprocessing,
    "signif_overlap",
    "Number of Significant Overlap for Category (permutation FDR < 0.05)",
    "Number of Significant Overlap",
    color_column = "colors_psychiatric",
    use_column_colors = TRUE  # Indicates that colors should be used based on the specified column
  )

# Plot 2: Visualize the number of significant phenotypes for each category
plot2 <-
  gene_overlap_barplot_summary(
    overlap_results_preprocessing,
    "signif_phenotypes",
    "Number of signif phenotypes (permutation FDR < 0.05)",
    "Number of Significant phenotypes",
    color_column = "colors_psychiatric",
    use_column_colors = TRUE
  )

# Plot 3: Visualize the number of significant GR (Gene-Related) lists for each category
plot3 <-
  gene_overlap_barplot_summary(
    overlap_results_preprocessing,
    "signif_papers",
    "Number signif GR lists (permutation FDR < 0.05)",
    "Number of Significant GR lists",
    color_column = "colors_psychiatric",
    use_column_colors = TRUE
  )

# Plot 4: Visualize the ratio of significant phenotypes to total phenotypes in each category
plot4 <-
  gene_overlap_barplot_summary(
    overlap_results_preprocessing,
    "ratio_signif_phenotypes",
    "Ratio signif phenotypes to phenotypes in category (permutation FDR < 0.05)",
    "(signif phenotypes)/(phenotypes in category)",
    color_column = "colors_psychiatric",
    use_column_colors = TRUE
  )

# Combine and display the plots
plot1 + plot2 + plot3 + plot4

################################################################################


# Prepare and visualize data with FDR and gene overlap count, coloring by FDR threshold with special handling for low counts
overlap_results_preprocessing %>% 
  mutate(fdr = map(data, ~ .x$original_data$df$fdr)) %>% 
  mutate(gene_overlap_count = map(data, ~ .x$original_data$df$gene_overlap_count)) %>% 
  select(category, fdr, FDR_threshold, gene_overlap_count) %>%
  unnest(c(fdr, gene_overlap_count)) %>% 
  mutate(FDR_threshold_colors = ifelse(FDR_threshold == "<0.001", "#F8766D",
                                       ifelse(FDR_threshold == ">0.05", "#619CFF", "#00BA38"))) %>% 
  mutate(FDR_threshold_colors = ifelse(gene_overlap_count < 3, "gray", FDR_threshold_colors)) %>% 
  mutate(category = as.factor(category)) %>%
  mutate(log10_fdr = -log10(fdr)) %>% 
  filter(log10_fdr < 25) %>%
  gene_overlap_manhattan_plot(color_column = "FDR_threshold_colors")

# Similar to the first block but skips filtering for log10_fdr and converting category to factor
overlap_results_preprocessing %>% 
  mutate(fdr = map(data, ~ .x$original_data$df$fdr)) %>% 
  mutate(gene_overlap_count = map(data, ~ .x$original_data$df$gene_overlap_count)) %>% 
  select(category, fdr, FDR_threshold, gene_overlap_count) %>%
  unnest(c(fdr, gene_overlap_count)) %>% 
  mutate(FDR_threshold_colors = ifelse(FDR_threshold == "<0.001", "#F8766D",
                                       ifelse(FDR_threshold == ">0.05", "#619CFF", "#00BA38"))) %>% 
  mutate(FDR_threshold_colors = ifelse(gene_overlap_count < 3, "gray", FDR_threshold_colors)) %>% 
  gene_overlap_manhattan_plot(color_column = "FDR_threshold_colors")

# Filter for significant overlaps (permutation FDR < 0.05), adjust colors for psychiatric association, and visualize
overlap_results_preprocessing %>% 
  filter(permutation_FDR < 0.05) %>% 
  mutate(gene_overlap_count = map(data, ~ .x$original_data$df$gene_overlap_count)) %>% 
  mutate(fdr = map(data, ~ .x$original_data$df$fdr)) %>% 
  unnest(c(gene_overlap_count, fdr))  %>% 
  mutate(psychiatric_associated = ifelse(colors_psychiatric == "black", "no", "yes")) %>% 
  mutate(colors_psychiatric = ifelse(gene_overlap_count < 3, "gray", colors_psychiatric)) %>% 
  mutate(colors_psychiatric = ifelse(colors_psychiatric == "black", "orange", colors_psychiatric)) %>% 
  select(category, fdr, psychiatric_associated, colors_psychiatric) %>%
  mutate(colors_psychiatric = as.factor(colors_psychiatric)) %>% 
  gene_overlap_manhattan_plot(data = ., color_column = "colors_psychiatric")


###############################################################################
# correlation analysis #
########################

# Calculate the number of unique genes associated with each category in the Biobank UK dataset
biobank_category_size_df <- map_dbl(categories_biobankuk, function(x) {
  # Filters phenotypes by category and counts the unique genes associated with them
  filter_phenotypes_by_category(
    genes_list = genes_phenotypes_PanUkBiobank,
    category_name = x,
    path_metafile = "data/phenotypes/panukbiobank-phenotype-category.csv"
  ) %>% 
    unname() %>% 
    unlist() %>% 
    unique() %>% 
    length()
}) %>% 
  enframe(name = "category", value = "n_genes")

# Analyze the overlap of genes between categories and prepare the data for visualization
biobank_category_overlap_info <- overlap_resutls_random_dependent_fdr0.01 %>%
  enframe(name = "category", value = "data") %>%
  # Calculate the number of unique overlapping genes for each category
  mutate(n_overlap_genes = map_dbl(data, ~ .x$significant_data$df$overlap_genes %>%
                                     strsplit(., ",") %>%
                                     unlist() %>%
                                     unique() %>%
                                     length())) %>%
  arrange(desc(n_overlap_genes)) %>%
  select(category, n_overlap_genes) %>%
  # Join the data with the category size information
  left_join(biobank_category_size_df, by = "category") %>%
  # Join with preprocessing results to get more data for plotting
  left_join(overlap_results_preprocessing[, -2], by = "category")

# Prepare the first plot showing the relationship between the number of genes and the number of significant phenotypes
p1 <- gene_overlap_regression_plot(
  data = biobank_category_overlap_info, 
  x_var = n_genes, 
  y_var = signif_phenotypes, 
  x_transform = TRUE, 
  y_transform = TRUE, 
  y_divisor = 1
)

# Prepare the second plot showing the relationship adjusted by the total number of phenotypes
p2 <- gene_overlap_regression_plot(
  data = biobank_category_overlap_info, 
  x_var = n_genes, 
  y_var = signif_phenotypes/number_phenotypes, 
  x_transform = TRUE, 
  y_transform = TRUE, 
  y_divisor = 1
)

# Combine and display the two plots side by side
p1 + p2
