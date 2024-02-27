# Define a function to create plots with categories ordered from largest to smallest
create_significance_plot <- function(data, y_value, title, y_lab) {
  # Ensure the colors are correctly mapped
  color_mapping <- setNames(data$colors_psychiatric, data$category)
  
  # Filter, mutate, and plot with categories ordered by y_value in descending order
  data %>%
    filter(permutation_FDR < 0.05) %>%
    mutate(category = reorder(category, get(y_value), FUN = function(x) -mean(x))) %>%
    ggplot(aes(x = category, y = get(y_value), fill = category)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = color_mapping, guide = FALSE) +
    theme_minimal() +
    labs(title = title,
         x = "Category",
         y = y_lab) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
}

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


# Use the function to create each plot
plot1 <- create_significance_plot(overlap_results_preprocessing, "signif_overlap", "Number of Significant Overlap for Category (permutation FDR < 0.05)", "Number of Significant Overlap")
plot2 <- create_significance_plot(overlap_results_preprocessing, "signif_phenotypes", "Number of signif phenotypes (permutation FDR < 0.05)", "Number of Significant phenotypes")
plot3 <- create_significance_plot(overlap_results_preprocessing, "signif_papers", "Number signif GR lists (permutation FDR < 0.05)", "Number of Significant GR lists")
plot4 <- create_significance_plot(overlap_results_preprocessing, "ratio_signif_phenotypes", "Ratio signif phenotypes to phenotypes in category (permutation FDR < 0.05)", "(signif phenotypes)/(phenotypes in category)")


plot1 + plot2 + plot3 + plot4

################################################################################
overlap_results_preprocessing %>% 
  mutate(fdr = map(data, ~ .x$original_data$df$fdr)) %>% 
  mutate(gene_overlap_count = map(data, ~ .x$original_data$df$gene_overlap_count)) %>% 
  select(category, fdr,  FDR_threshold, gene_overlap_count) %>%
  unnest(c(fdr, gene_overlap_count)) %>% 
  mutate(FDR_threshold_colors = ifelse(FDR_threshold == "<0.001", "#F8766D",
                                       ifelse(FDR_threshold == ">0.05", "#619CFF", "#00BA38"))) %>% 
  mutate(FDR_threshold_colors = ifelse(gene_overlap_count < 3, "gray", FDR_threshold_colors)) %>% 
  # Convert category to a factor for better plotting
  mutate(category = as.factor(category)) %>%
  mutate(log10_fdr = -log10(fdr)) %>% 
  filter(log10_fdr < 25) %>%
  gene_overlap_manhattan_plot(color_column = "FDR_threshold_colors")
  



overlap_results_preprocessing %>% 
  mutate(fdr = map(data, ~ .x$original_data$df$fdr)) %>% 
  mutate(gene_overlap_count = map(data, ~ .x$original_data$df$gene_overlap_count)) %>% 
  # Select only the necessary columns and unnest the chi2 values
  select(category, fdr, FDR_threshold, gene_overlap_count) %>%
  unnest(c(fdr, gene_overlap_count)) %>% 
  mutate(FDR_threshold_colors = ifelse(FDR_threshold == "<0.001", "#F8766D",
                                       ifelse(FDR_threshold == ">0.05", "#619CFF", "#00BA38"))) %>% 
  mutate(FDR_threshold_colors = ifelse(gene_overlap_count < 3, "gray", FDR_threshold_colors)) %>% 
  gene_overlap_manhattan_plot()
  
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
  