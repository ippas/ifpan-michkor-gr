library(ggplot2)
library(reshape2)

# Sample data
data <- data.frame(
  Gene = c("APOB", "APOC3", "LDLR", "ABCA1", "LIPC", "KLC3", "SH2B3", "HECTD4"),
  Phenotype = c(1, 1, 1, 0, 1, 0, 0, 1),
  Expression = c(1, 0, 1, 1, 0, 1, 0, 0),
  Localization = c(0, 1, 0, 0, 1, 1, 1, 0),
  Pathway = c(1, 0, 1, 1, 0, 1, 0, 1)
)

data

# Melting the data
melted_data <- melt(data, id.vars = "Gene")

melted_data

# Create the heatmap
ggplot(melted_data, aes(x = variable, y = Gene, fill = factor(value))) +
  geom_tile(color = "white") + # Add tiles with white borders
  scale_fill_manual(values = c("0" = "grey", "1" = "red")) + # Set custom colors
  theme_minimal() + # Use a minimal theme
  theme(axis.text.x = element_text(angle = 90, hjust = 1), # Rotate x-axis labels for better visibility
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank()) # Remove minor grid lines


# Create the heatmap with adjusted styles
ggplot(melted_data, aes(x = variable, y = Gene, fill = factor(value))) +
  geom_tile(color = "black", size = 1.5) + # Black borders with increased size for visibility
  scale_fill_manual(values = c("0" = "grey", "1" = "red")) + # Set custom colors
  theme_minimal() + # Use a minimal theme
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), # Rotate x-axis labels
        axis.title.x = element_blank(), # Hide x-axis title
        axis.title.y = element_blank(), # Hide y-axis title
        panel.spacing = unit(2, "lines"), # Increase spacing between cells
        panel.background = element_rect(fill = "white", colour = "black"), # White background with black outline
        plot.background = element_rect(fill = "white", colour = NA), # Transparent plot background
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank(), # Remove minor grid lines
        strip.background = element_rect(fill = "white", colour = "black")) # Adjust the background of facet labels

# Create the heatmap with adjusted styles
ggplot(melted_data, aes(x = variable, y = Gene, fill = factor(value))) +
  geom_tile(color = "black", size = 1.5, width = 0.9, height = 0.9) + # Adjust width and height for spacing
  scale_fill_manual(values = c("0" = "grey", "1" = "red")) + # Set custom colors
  theme_minimal() + # Use a minimal theme
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), # Rotate x-axis labels
        axis.title.x = element_blank(), # Hide x-axis title
        axis.title.y = element_blank(), # Hide y-axis title
        panel.spacing = unit(4, "lines"), # Increase spacing between cells
        panel.background = element_rect(fill = "white", colour = "black"), # White background with black outline
        plot.background = element_rect(fill = "white", colour = NA), # Transparent plot background
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank(), # Remove minor grid lines
        strip.background = element_rect(fill = "white", colour = "black")) # Adjust the background of facet labels


create_heatmap <- function(data, panel_spacing_size = 4, border_size = 1.5) {
  # Assuming data is already in the correct format (melted if necessary)
  
  # Create the heatmap
  ggplot(data, aes(x = variable, y = Gene, fill = factor(value))) +
    geom_tile(color = "black", size = border_size, width = 0.9, height = 0.9) + # Tiles with custom border
    scale_fill_manual(values = c("0" = "grey", "1" = "red")) + # Custom colors
    theme_minimal() + # Minimal theme
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), # Rotate x-axis labels
          axis.title.x = element_blank(), # Hide x-axis title
          axis.title.y = element_blank(), # Hide y-axis title
          pan  geom_tile(data = to_plot, aes(x = variable, y = cluster, fill = factor(value)), color = "black", size = 1, width = 0.7, height = 0.7, na.rm = TRUE) +
            scale_fill_manual(values = c("0" = "#e2e2e2", "1" = "#ff964f", "2" = "#db5856", "3" = "#592D1D")) + # Custom colorsel.spacing = unit(panel_spacing_size, "lines"), # Custom spacing between cells
          panel.background = element_rect(fill = "white", colour = "black"), # White background with black outline
          plot.background = element_rect(fill = "white", colour = NA), # Transparent plot background
          panel.grid.major = element_blank(), # Remove major grid lines
          panel.grid.minor = element_blank(), # Remove minor grid lines
          strip.background = element_rect(fill = "white", colour = "black")) # Adjust the background of facet labels
}

create_heatmap <- function(data, panel_spacing_size = 4, border_size = 1.5, tile_size = 0.9) {
  # Create the heatmap
  ggplot(data, aes(x = variable, y = Gene, fill = factor(value))) +
    geom_tile(color = "black", size = border_size, width = tile_size, height = tile_size) + # Tiles with custom sizes
    scale_fill_manual(values = c("0" = "grey", "1" = "red")) + # Custom colors
    theme_minimal() + # Minimal theme
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), # Rotate x-axis labels
          axis.title.x = element_blank(), # Hide x-axis title
          axis.title.y = element_blank(), # Hide y-axis title
          panel.spacing = unit(panel_spacing_size, "lines"), # Custom spacing between cells
          panel.background = element_rect(fill = "white", colour = "black"), # White background with black outline
          plot.background = element_rect(fill = "white", colour = NA), # Transparent plot background
          panel.grid.major = element_blank(), # Remove major grid lines
          panel.grid.minor = element_blank(), # Remove minor grid lines
          strip.background = element_rect(fill = "white", colour = "black")) # Adjust the background of facet labels
}

create_heatmap <- function(data, panel_spacing_size = 6, border_size = 1.5, tile_size = 0.7) {
  # Create the heatmap
  plot <- ggplot(data, aes(x = variable, y = Gene, fill = factor(value))) +
    geom_tile(color = "black", size = border_size, width = tile_size, height = tile_size) + # Tiles with custom sizes
    scale_fill_manual(values = c("0" = "grey", "1" = "red")) + # Custom colors
    theme_minimal() + # Minimal theme
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), # Rotate x-axis labels
          axis.title.x = element_blank(), # Hide x-axis title
          axis.title.y = element_blank(), # Hide y-axis title
          panel.spacing = unit(panel_spacing_size, "lines"), # Custom spacing between cells
          plot.background = element_rect(fill = "white", colour = NA), # Transparent plot background
          panel.grid.major = element_blank(), # Remove major grid lines
          panel.grid.minor = element_blank(), # Remove minor grid lines
          strip.background = element_rect(fill = "white", colour = "black"), # Adjust the background of facet labels
          panel.border = element_blank()) # Remove default panel border
  
  # Add a vertical line on the left side only
  plot + geom_vline(xintercept = 0.5, color = "black", size = 1)
}
create_heatmap <- function(data, panel_spacing_size = 6, border_size = 1.5, tile_size = 0.7) {
  # Adjusting the factor levels to add space between specific columns
  data$variable <- factor(data$variable, levels = c("Phenotype", "blank", "Expression", "Localization", "Pathway"))
  
  # Create the heatmap
  plot <- ggplot(data, aes(x = variable, y = Gene, fill = factor(value))) +
    geom_tile(color = "black", size = border_size, width = tile_size, height = tile_size, na.rm = TRUE) + # Tiles with custom sizes
    scale_fill_manual(values = c("0" = "grey", "1" = "red")) + # Custom colors
    theme_minimal() + # Minimal theme
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), # Rotate x-axis labels
          axis.title.x = element_blank(), # Hide x-axis title
          axis.title.y = element_blank(), # Hide y-axis title
          panel.spacing = unit(panel_spacing_size, "lines"), # Custom spacing between cells
          plot.background = element_rect(fill = "white", colour = NA), # Transparent plot background
          panel.grid.major = element_blank(), # Remove major grid lines
          panel.grid.minor = element_blank(), # Remove minor grid lines
          strip.background = element_rect(fill = "white", colour = "black"), # Adjust the background of facet labels
          panel.border = element_blank()) # Remove default panel border
  
  # Add a vertical line on the left side only
  plot + geom_vline(xintercept = 0.5, color = "black", size = 1)
}
create_heatmap(data = melted_data, border_size = 1, tile_size = 0.7)

create_heatmap <- function(data, x_axis, y_axis, panel_spacing_size = 6, border_size = 1.5, tile_size = 0.7,
                           x_axis_text_size = 12, y_axis_text_size = 12) {
  # Sort y-axis elements in reverse alphabetical order and adjust as a factor
  data <- data %>%
    arrange(desc(!!sym(y_axis))) %>%
    mutate(!!y_axis := factor(!!sym(y_axis), levels = unique(!!sym(y_axis))))
  
  # Create the heatmap
  plot <- ggplot(data, aes(x = get(x_axis), y = get(y_axis), fill = factor(value))) +
    geom_tile(color = "black", size = border_size, width = tile_size, height = tile_size, na.rm = TRUE) + # Tiles with custom sizes
    scale_fill_manual(values = c("0" = "#e2e2e2", "1" = "#300000", "2" = "#ff4122", "3" = "#c61a09")) + # Custom colors
    # scale_fill_manual(values = c("0" = "grey", "1" = "orange", "2" = "red")) + # Custom colors
    theme_minimal() + # Minimal theme
    theme(
      axis.text.x.top = element_text(angle = 90, hjust = 1, size = x_axis_text_size, color = "black", vjust = 0.5), # X-axis labels at the top
      axis.ticks.x.top = element_line(color = "black"), # Ensure ticks are at top if needed
      axis.text.x = element_blank(), # Hide the default x-axis text at bottom
      axis.text.y = element_text(size = y_axis_text_size), # Y-axis text size
      axis.title.x.top = element_blank(), # Ensure no title at the top x-axis
      axis.title.x = element_blank(), # Hide the bottom x-axis title if it exists
      axis.title.y = element_blank(), # Hide y-axis title
      axis.ticks.length = unit(0, "points"), # Remove tick marks
      panel.spacing = unit(panel_spacing_size, "lines"), # Custom spacing between cells
      plot.background = element_rect(fill = "white", colour = NA), # Transparent plot background
      panel.grid.major = element_blank(), # Remove major grid lines
      panel.grid.minor = element_blank(), # Remove minor grid lines
      strip.background = element_rect(fill = "white", colour = "black"), # Adjust the background of facet labels
      panel.border = element_blank(), # Remove default panel border
      legend.position = "bottom", # Move legend to bottom
      legend.box = "horizontal", # Align legend items horizontally
      legend.title = element_blank() # Remove legend title if desired
    ) +
    scale_x_discrete(position = "top")  # Move the x-axis labels to the top
  
  # Add a vertical line on the left side only
  plot + geom_vline(xintercept = 0.5, color = "black", size = 1)
}
