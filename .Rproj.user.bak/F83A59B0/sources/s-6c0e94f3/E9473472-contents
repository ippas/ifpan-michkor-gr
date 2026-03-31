create_gene_pie_chart <- function(gene_list, threshold_percent = 1) {
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
    other_row <- data.frame(frequency = paste0("other (> ", data[nrow(data), 1], ")"), count = other_count, percent = other_percent)
    data <- bind_rows(data, other_row)
  }
  
  print(data)
  # Plot the pie chart using ggplot2 based on the percentage column
  p <- ggplot(data, aes(x = "", y = percent, fill = factor(frequency))) +
    geom_bar(width = 1, stat = "identity", color = "white") +
    coord_polar("y", start = 0) +
    theme_void() +
    labs(fill = "Frequency Count") +
    geom_text(aes(label = paste(round(percent, 1), "%")), position = position_stack(vjust = 0.5)) +
    # theme(legend.position = "bottom") +  # Move legend to bottom
    guides(fill = guide_legend(title = "Gene Frequency", ncol = 1))  
  
  return(p)
}

# 
# create_gene_pie_chart(gene_list = papers_gene_list, threshold_percent = 96) 
# 
# papers_gene_list[grep("_down", names(papers_gene_list))] %>% length()
# 
# 
# 
# create_gene_pie_chart(gene_list = papers_gene_list[grep("_down", names(papers_gene_list))], threshold_percent = 96) 
# 
# 
# create_gene_pie_chart(gene_list = papers_gene_list[grep("_up", names(papers_gene_list))], threshold_percent = 96) 
