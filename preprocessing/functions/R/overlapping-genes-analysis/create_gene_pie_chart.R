create_gene_pie_chart <- function(gene_list, threshold_percent = 99, palette = "viridis") {
  # Convert the list to a vector and create a table of gene occurrences
  gene_frequencies <- table(unlist(gene_list))
  
  # Calculate how often each frequency occurs
  frequency_of_frequencies <- table(gene_frequencies)
  
  # Convert this table to a data frame
  data <- as.data.frame(frequency_of_frequencies)
  
  # Set column names to snake_case
  colnames(data) <- c("frequency", "count")
  
  # Calculate the percentage of each frequency count
  data$percent <- (data$count / sum(data$count)) * 100
  
  # Calculate the cumulative sum of percentages
  data$cumulative_percent <- cumsum(data$percent)
  
  # Determine the threshold count based on the cumulative sum
  threshold_count <- sum(data$count[data$cumulative_percent < threshold_percent])
  
  # Group frequencies below threshold_count into "Other"
  if (threshold_percent < 100) {
    other_count <- sum(data$count[data$cumulative_percent >= threshold_percent])
    other_percent <- sum(data$percent[data$cumulative_percent >= threshold_percent])
    data <- data %>% filter(cumulative_percent < threshold_percent)
    # other_row <- data.frame(frequency = paste0("other (> ", data[nrow(data), 1], ")"), count = other_count, percent = other_percent)
    other_row <- data.frame(frequency = paste0("more"), count = other_count, percent = other_percent)
    data <- bind_rows(data, other_row)
  }
  
  # Print the data for verification
  print(data)
  
  # Plot the pie chart using ggplot2 based on the percentage column
  p <- ggplot(data, aes(x = "", y = percent, fill = factor(frequency))) +
    geom_bar(width = 1, stat = "identity", color = "white") +
    coord_polar("y", start = 0) +
    theme_void() +
    labs(fill = "Frequency Count") +
    geom_text(aes(label = paste(round(percent, 1), "%")), position = position_stack(vjust = 0.5)) +
    guides(fill = guide_legend(title = "Gene Frequency", ncol = 1)) +
    scale_fill_manual(values = palette) # Use the specified palette
  
  return(p)
}


generate_palette <- function(n, name) {
  # Ensure n is within the allowed range for Pastel1
  if (n > 9) {
    stop("n must be less than or equal to 9 for 'Pastel1'")
  }
  
  # Generate palette using brewer.pal
  colors <- brewer.pal(n, name)
  
  # Create named vector for colors
  names(colors) <- as.character(seq_len(n))
  
  # Add "more" color
  colors <- c(colors)
  
  # Change the name of the last element
  names(colors)[length(colors)] <- "more"
  
  return(colors)
}
