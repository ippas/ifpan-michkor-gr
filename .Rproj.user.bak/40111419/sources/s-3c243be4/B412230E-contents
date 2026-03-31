gene_overlap_barplot_summary <- function(data, y_value, title = "Significance Plot", y_lab = "Value", color_column = NULL, use_column_colors = FALSE, fdr_threshold = 0.05, text_angle = 45) {

  # Function: gene_overlap_barplot_summary
  # Description: Creates a bar plot to summarize gene overlap data, highlighting significance based on a specified y-value.
  #              The function allows for dynamic color mapping based on a specified column and can filter data based on an FDR threshold.
  #              Users can choose to use specified colors in the color_column or generate a color palette.
  #
  # Parameters:
  #   data - A dataframe containing the gene overlap data. Must include columns specified by y_value and optionally color_column.
  #   y_value - The name of the column in 'data' to be used for the y-axis values in the bar plot.
  #   title - A string specifying the title of the plot. Defaults to "Significance Plot".
  #   y_lab - The label for the y-axis. Defaults to "Value".
  #   color_column - Optional. The name of the column in 'data' used for coloring the bars. If NULL, bars are colored gray by default.
  #   use_column_colors - Boolean. If TRUE, uses the colors specified in color_column directly. If FALSE, generates a color palette.
  #   fdr_threshold - A numeric value specifying the FDR threshold for filtering the data. Defaults to 0.05.
  #   text_angle - The angle for the x-axis text labels. Defaults to 45 degrees.
  #
  # Returns:
  #   A ggplot object representing the bar plot of gene overlap data, with bars colored based on the specified color_column or a generated palette.
  
  # Validate the presence of the y_value column in the data
  if (!y_value %in% names(data)) {
    stop("The specified y_value column does not exist in the data frame.")
  }
  
  # Optionally filter the data based on the FDR threshold
  if ("permutation_FDR" %in% names(data)) {
    data <- data %>% filter(permutation_FDR < fdr_threshold)
  }
  
  # Prepare the plot with basic aesthetics
  plot <- data %>%
    mutate(category = reorder(category, get(y_value), FUN = function(x) -mean(x))) %>%
    ggplot(aes(x = category, y = get(y_value))) +
    geom_bar(stat = "identity") +  # Use actual values for bar heights
    theme_minimal() +  # Apply a minimal theme for a clean appearance
    labs(title = title, x = "Category", y = y_lab) +  # Set custom plot title and axis labels
    theme(axis.text.x = element_text(angle = text_angle, hjust = 1))  # Rotate x-axis labels for better readability
  
  # Conditionally apply colors based on the color_column
  if (!is.null(color_column) && color_column %in% names(data)) {
    if (use_column_colors) {
      # Directly use the colors specified in the color_column
      plot <- plot + aes(fill = !!sym(color_column)) + scale_fill_identity()
    } else {
      # Generate a color palette if not using specified colors
      unique_colors <- length(unique(data[[color_column]]))
      color_values <- RColorBrewer::brewer.pal(min(unique_colors, 8), "Set3")
      plot <- plot + aes(fill = factor(!!sym(color_column))) + scale_fill_manual(values = color_values)
    }
  } else {
    # Apply a default gray color if no color_column is specified
    plot <- plot + scale_fill_manual(values = "gray")
  }
  
  return(plot)
}
