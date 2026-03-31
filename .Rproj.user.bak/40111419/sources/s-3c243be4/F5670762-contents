# gene_overlap_manhattan_plot <- function(data, color_column = NULL, title = "Gene Overlap Manhattan Plot", y_intercepts = c(2, 10, 20)) {
#   # Function: gene_overlap_manhattan_plot
#   # Description: Creates a scatter plot resembling a Manhattan plot to visualize gene overlap data. 
#   # This function plots -log10 transformed FDR values across different categories and allows for optional coloring based on a specified column. 
#   # Horizontal lines can be added at specified y-values to indicate significance levels.
#   #
#   # Parameters:
#   #   data - A dataframe containing gene overlap data. Must include 'category' and 'fdr' columns.
#   #   color_column - Optional. A string specifying the column name to be used for coloring the points. If NULL, no color differentiation is applied.
#   #   title - The title of the plot. Defaults to "Gene Overlap Manhattan Plot".
#   #   y_intercepts - A numeric vector indicating y-values at which horizontal lines should be drawn to denote different significance levels. Defaults to c(2, 10, 20).
#   #
#   # Returns:
#   #   A ggplot object that visually represents the gene overlap data, suitable for exploratory data analysis and presentation.
#   
#   # Ensure 'category' is treated as a factor for proper plotting
#   data <- data %>%
#     mutate(category = as.factor(category),
#            log10_fdr = -log10(fdr)) # Calculate the -log10 of FDR values for a more interpretable visualization scale
#   
#   # Initialize the ggplot with basic aesthetics
#   p <- ggplot(data, aes(x = category, y = log10_fdr)) +
#     geom_jitter() +  # Add jitter to points to reduce overlap and improve visibility
#     theme_minimal() +  # Use a minimal theme for a clean plot appearance
#     labs(title = title, x = "Category", y = "-log10(FDR)") +  # Set plot title and axis labels
#     theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability
#   
#   # Conditionally add color mapping if a valid color_column is provided
#   if (!is.null(color_column) && color_column %in% names(data)) {
#     p <- p + aes_string(color = color_column) + scale_color_identity()  # Dynamically map the color aesthetic to the specified column
#   }
#   
#   # Add horizontal lines at specified y-intercepts to denote significance levels
#   p <- p + map(y_intercepts, ~geom_hline(yintercept = .x, linetype = "dashed", color = "black", size = 1))
#   
#   return(p)
# }

gene_overlap_manhattan_plot <- function(data, color_column = NULL, title = "Gene Overlap Manhattan Plot", y_intercepts = c(2, 10, 20)) {
  # Ensure 'category' is treated as a factor for proper plotting
  data <- data %>%
    mutate(category = as.factor(category),
           log10_fdr = -log10(fdr)) # Calculate the -log10 of FDR values for a more interpretable visualization scale
  
  # Convert color_column to symbol if not NULL
  color_col <- ensym(color_column)
  
  # Initialize the ggplot
  p <- ggplot(data, aes(x = category, y = log10_fdr)) +
    geom_jitter() +  # Add jitter to points to reduce overlap and improve visibility
    theme_minimal() +  # Use a minimal theme for a clean plot appearance
    labs(title = title, x = "Category", y = "-log10(FDR)") +  # Set plot title and axis labels
    theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for better readability
          legend.position = "bottom")  # Position the legend at the bottom
  
  # Conditionally add color mapping
  if (!is.null(color_column)) {
    p <- p + aes(color = !!color_col)
  }
  
  # Add horizontal lines at specified y-intercepts to denote significance levels
  p <- p + map(y_intercepts, ~geom_hline(yintercept = .x, linetype = "dashed", color = "black", size = 1))
  
  return(p)
}
