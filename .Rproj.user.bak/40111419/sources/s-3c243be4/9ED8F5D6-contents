library(umap)
library(ggplot2)

# Function to apply UMAP with combinations of n_neighbors and min_dist, PCA, and generate plots
# Now with labels as an optional parameter
run_dimensionality_reduction_combinations_with_plots <- function(data_matrix, n_neighbors_vector, min_dist_vector, labels = NULL) {
  results <- list()
  
  # Perform PCA for comparison and plot
  pca_result <- prcomp(data_matrix, center = TRUE, scale. = TRUE)
  pca_data <- as.data.frame(pca_result$x[, 1:2])
  pca_data$labels <- if (!is.null(labels)) labels else rep("All", nrow(pca_data))
  pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = labels)) +
    geom_point() +
    labs(title = "PCA", x = "PC1", y = "PC2") +
    theme_minimal() +
    scale_color_discrete(name = if (is.null(labels)) NULL else "Labels")
  results$pca <- list(embedding = pca_result$x[, 1:2], plot = pca_plot)
  
  # Generate all combinations of n_neighbors and min_dist
  param_combinations <- expand.grid(n_neighbors = n_neighbors_vector, min_dist = min_dist_vector)
  
  # Iterate over each combination of parameters
  for (i in 1:nrow(param_combinations)) {
    params <- param_combinations[i, ]
    set.seed(42) # For reproducibility
    umap_result <- umap(data_matrix, n_neighbors = params$n_neighbors, min_dist = params$min_dist)
    
    # Prepare data for plotting
    umap_data <- as.data.frame(umap_result$layout)
    colnames(umap_data) <- c("UMAP1", "UMAP2")
    umap_data$labels <- if (!is.null(labels)) labels else rep("All", nrow(umap_data))
    umap_plot <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = labels)) +
      geom_point() +
      labs(title = paste("UMAP with n_neighbors =", params$n_neighbors, "& min_dist =", params$min_dist), x = "UMAP1", y = "UMAP2") +
      theme_minimal() +
      scale_color_discrete(name = if (is.null(labels)) NULL else "Labels")
    
    # Store the UMAP result and plot
    results[[paste("umap_n", params$n_neighbors, "_mdist", params$min_dist, sep = "")]] <- list(embedding = umap_result$layout, plot = umap_plot)
  }
  
  return(results)
}

# Define vectors of parameters to explore
n_neighbors_vector <- c(5, 15, 30, 50, 60, 70, 80, 100)
min_dist_vector <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)

# Example usage without labels
# Assuming genes_phenotypes_PanUkBiobank_matrix is your data matrix
results2 <- run_dimensionality_reduction_combinations_with_plots(genes_phenotypes_PanUkBiobank_matrix, n_neighbors_vector, min_dist_vector)
# Access and display a specific plot
# For example, to display the plot for UMAP with n_neighbors = 15 and min_dist = 0.1:
print(results$umap_n50_mdist0.5$plot)



# Assuming the rest of the function is the same and focusing on the plotting part
combine_umap_plots <- function(results) {
  # Initialize a list to store the ggplot objects
  plot_list <- list()
  
  # Iterate through the results to add UMAP plots to the list
  for (result_name in names(results)) {
    if (result_name != "pca") { # Optionally exclude PCA from the combined plot
      plot_list[[result_name]] <- results[[result_name]]$plot
    }
  }
  
  # Use wrap_plots() to combine all plots in the list
  # The layout can be adjusted by specifying the number of rows (nrow) and columns (ncol)
  combined_plot <- wrap_plots(plotlist = plot_list, ncol = 4) # Adjust ncol and nrow as needed
  
  return(combined_plot)
}


combine_umap_plots(results)

results$umap_n50_mdist0.5$embedding %>% 
  as.data.frame() %>% 
  set_colnames(c("UMAP1", "UMAP2")) %>% 
  # filter(UMAP1 > 25) 
  filter(UMAP2 < 0) 
  ggplot(aes(x = UMAP1, y = UMAP2)) +
  geom_point() + 
  theme_minimal()
