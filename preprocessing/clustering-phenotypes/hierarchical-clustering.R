calculate_distance_matrices <- function(data_matrix, minkowski_p = 3) {
  # Load necessary library, if not already loaded
  if (!requireNamespace("stats", quietly = TRUE)) install.packages("stats")
  library(stats)
  
  # Define the list of methods to use, excluding Minkowski for now
  methods <- c("euclidean", "maximum", "manhattan", "canberra", "binary")
  
  # Initialize an empty list to store distance matrices
  distance_matrices <- list()
  
  # Loop through each method and calculate the distance matrix
  for (method in methods) {
    print(method)
    dist_matrix <- dist(data_matrix, method = method)
    # Add the distance matrix to the list with the method name as the key
    distance_matrices[[method]] <- dist_matrix
  }
  
  print("minkowski")
  # Calculate the Minkowski distance matrix with specified power and add it to the list
  minkowski_dist_matrix <- dist(data_matrix, method = "minkowski", p = minkowski_p)
  distance_matrices[["minkowski"]] <- minkowski_dist_matrix
  
  # Return the list of distance matrices
  return(distance_matrices)
}

# Function to run hclust with all methods
run_hclust_all_methods <- function(dist_matrix) {
  # Define the methods to be used in hclust
  methods <- c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")
  
  # Initialize an empty list to store hclust results
  hclust_results <- list()
  
  # Loop through each method and perform hierarchical clustering
  for (method in methods) {
    # Perform hierarchical clustering using the current method
    hclust_result <- hclust(dist_matrix, method = method)
    
    # Store the result in the list with method name as the key
    hclust_results[[method]] <- hclust_result
  }
  
  # Return the list of hclust results
  return(hclust_results)
}


genes_phenotypes_PanUkBiobank_matrix -> genes_phenotypes_PanUkBiobank_matrix2

rownames(genes_phenotypes_PanUkBiobank_matrix2) <- gsub("biobankuk-", "", rownames(genes_phenotypes_PanUkBiobank_matrix2))

phenotypes_dist_list <- calculate_distance_matrices(genes_phenotypes_PanUkBiobank_matrix2)


# Plot the dendrogram
plot(hc, hang = -1)
hc <- hclust(phenotypes_dist_list$binary, method = "mcquitty")

phenotypes_dist_list %>% 
  lapply(., run_hclust_all_methods) ->  hclust_results_list


tmp_list %>%
  # Apply cutree and table to each element in the list
  map(~ cutree(.x, k = 200) %>% table() %>% as.vector()) %>%
  # Convert the list to a data frame
  map_dfr(~ as.data.frame(t(.)), .id = "Method") %>%
  as.data.frame() %>%
  gather(., key = "colname", value = "number", -Method) %>%  
  select(-c(colname)) %>% 
  mutate(log2_number = log2(number)) %>% 
  ggplot(aes(x = Method, y = log2_number)) + 
  geom_boxplot(color = "gray") +
  geom_jitter() +
  theme_minimal()

# Updated function with 'k' as an argument
process_hclust_results <- function(hclust_list, phenotype_name, k) {
  hclust_list %>%
    map(~ cutree(.x, k = k) %>% table() %>% as.vector()) %>%
    map_dfr(~ as.data.frame(t(.)), .id = "Method") %>%
    mutate(Phenotype = phenotype_name) %>%
    gather(key = "colname", value = "number", -Method, -Phenotype) %>%
    select(-colname) %>%
    mutate(log2_number = log2(number))
}

k_value <- 200

# Apply the function to each element in hclust_results_list and combine into one data frame
all_phenotypes_df <- map2_df(hclust_results_list, names(hclust_results_list), process_hclust_results, k_value)

# Now plot with facet_wrap
ggplot(all_phenotypes_df, aes(x = Method, y = log2_number)) + 
  geom_jitter() +
  theme_minimal() +
  geom_boxplot(color = "orange", fill = NA) +
  facet_wrap(~ Phenotype, scales = "free_y") # Add facet wrap for each phenotype  


################################################################################
# Define the list of k values
k_values <- c(10, 50, 100, 200, 500, 1000)

# Function to process and plot results for each k value
process_and_plot <- function(k_value) {
  # Apply the function to each element in hclust_results_list and combine into one data frame
  all_phenotypes_df <- map2_df(hclust_results_list, names(hclust_results_list), ~process_hclust_results(.x, .y, k_value))
  
  # Plot with facet_wrap
  p <- ggplot(all_phenotypes_df, aes(x = Method, y = log2_number)) + 
    geom_jitter() +
    theme_minimal() +
    geom_boxplot(color = "orange", fill = NA) +
    facet_wrap(~ Phenotype, scales = "free_y") + # Add facet wrap for each phenotype
    labs(title = paste("Results for k =", k_value)) # Add a title to indicate the k value
  
  print(p) # Print the plot
  
  return(p)
}

# Apply the plotting function for each k value
lapply(k_values, process_and_plot) -> sub_group_hclust_plot

################################################################################

# Open PDF device
pdf("results/figures-phenotypes/phenotypes-hclust.pdf", width = 1000, height = 50)

# Plotting (replace with your actual plot code)
plot(hc, hang = -1) # Make sure 'hc' is defined and contains your hierarchical clustering
rect.hclust(hc, k = 10, border = 2:10)
# Close PDF device
dev.off()


sub_grp <- cutree(hc, k = 10)

table(sub_grp)



hc <- hclust(dist_matrix, method = "complete")