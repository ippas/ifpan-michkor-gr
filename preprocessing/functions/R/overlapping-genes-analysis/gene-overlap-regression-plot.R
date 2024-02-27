gene_overlap_regression_plot <- function(data, x_var, y_var, x_transform = TRUE, y_transform = TRUE, y_divisor = 1) {
  
  # gene_overlap_regression_plot function
  # This function creates a regression plot for two variables from a given dataset, optionally applying log2 transformations
  # and adjusting the y variable by a divisor. It also calculates and displays the correlation coefficient on the plot.
  # 
  # Args:
  #   data: A data frame containing the data to be plotted.
  #   x_var: Unquoted name of the x variable or an expression involving column names.
  #   y_var: Unquoted name of the y variable or an expression involving column names.
  #   x_transform: Boolean indicating whether to apply a log2 transformation to the x variable.
  #   y_transform: Boolean indicating whether to apply a log2 transformation to the y variable.
  #   y_divisor: Numeric value to divide the y variable, default is 1 (no division).
  # 
  # Returns:
  #   A ggplot object representing the regression plot with annotated correlation coefficient.
  # Capture the expressions for x_var and y_var, and also convert them to strings for axis labels.
  
  x_var_expr <- enquo(x_var)
  y_var_expr <- enquo(y_var)
  x_var_label <- deparse(substitute(x_var)) # Converts the x_var expression into a string for labeling.
  y_var_label <- deparse(substitute(y_var)) # Converts the y_var expression into a string for labeling.
  
  # Evaluate the captured expressions within the context of the data frame, adjusting y_var by the divisor.
  data <- data %>%
    mutate(x_var_eval = eval_tidy(x_var_expr, .), # Evaluates x_var expression in the data frame context.
           y_var_eval = eval_tidy(y_var_expr, .) / y_divisor) # Evaluates y_var expression and applies divisor.
  
  # Apply log2 transformation if indicated by parameters and update axis labels to reflect transformations.
  if (x_transform) {
    data$x_var_eval <- log2(data$x_var_eval)
    x_var_label <- paste0("log2(", x_var_label, ")") # Updates the label to indicate log2 transformation.
  }
  
  if (y_transform) {
    data$y_var_eval <- log2(data$y_var_eval)
    y_var_label <- paste0("log2(", y_var_label, ")") # Updates the label to indicate log2 transformation.
  }
  
  # Calculate the correlation coefficient between the evaluated x and y variables.
  correlation <- cor(data$x_var_eval, data$y_var_eval)
  
  # Create the ggplot object, setting up aesthetics, points, regression line, and annotations.
  p <- ggplot(data, aes(x = x_var_eval, y = y_var_eval)) +
    geom_point() + # Adds points to the plot.
    geom_smooth(method = "lm", se = TRUE, color = "black") + # Adds a linear regression line with confidence interval.
    annotate("text", x = Inf, y = Inf, label = paste("Correlation:", round(correlation, 2)), 
             hjust = 1.1, vjust = 2, size = 5, color = "red") + # Annotates plot with correlation coefficient.
    labs(x = x_var_label, y = y_var_label) + # Sets custom axis labels based on input expressions.
    theme_minimal() # Uses a minimal theme for aesthetics.
  
  return(p) # Returns the ggplot object.
}

