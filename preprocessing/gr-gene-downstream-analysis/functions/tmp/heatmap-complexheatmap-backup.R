library(ComplexHeatmap)
library(circlize)

if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  install.packages("RColorBrewer")
}
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(grid)

draw_custom_heatmap <- function(data_list, 
                                data_type,
                                fdr_threshold = 0.05,
                                color_rect = "green",
                                lwd_rect = 2,
                                alpha_rect = 1,
                                apply_filling = TRUE,
                                color_filling = "green",
                                size_filling = 1,
                                alpha_filling = 1,
                                pch_filling = 16,
                                palette = c(
                                  "pastel_blue"       = "#69A6D1",
                                  "pastel_light_blue" = "#94DFFF",
                                  "white"        = "#ffffcc",
                                  "pastel_orange"= "#fed976",
                                  "pastel_red"   = "#fc4e2a"
                                ),
                                col_mapping_vector=NULL, 
                                row_mapping_vector=NULL,
                                col_significant=F
) {
  
  # Define the function to filter matrices
  # This internal function takes a list of matrices and filters each matrix based on the specified rows and columns.
  filter_matrices <- function(data_list, rows_to_filter, cols_to_filter) {
    # Initialize an empty list to store the filtered matrices
    filtered_list_significant <- list()
    
    # Loop through each item in the matrix list by their names
    for (name in names(data_list)) {
      # Extract the matrix using its name
      matrix <- data_list[[name]]
      # Apply the filter using the provided indices for rows and columns
      filtered_matrix <- matrix[rows_to_filter, cols_to_filter]
      # Add the filtered matrix to the list with the same name
      filtered_list_significant[[name]] <- filtered_matrix
    }
    
    # Return the list of filtered matrices
    return(filtered_list_significant)
  }
  
  data_list[[data_type]]$list -> data
  
  data_list[[data_type]]$cols -> cols_to_filter
  data_list[[data_type]]$rows -> rows_to_filter
  
  if(col_significant){
    filter_matrices(data = data,
                    rows_to_filter = rows_to_filter,
                    cols_to_filter = cols_to_filter) -> data
  }
  
  replace_names_in_data(data = data, 
                        col_mapping_vector = col_mapping_vector,
                        row_mapping_vector = row_mapping_vector) -> data
  
  chi2_matrix <- data$chi2_value_matrix
  chi2_matrix <- log2(chi2_matrix + 1)
  
  number_overlap_matrix <- data$number_overlap_matrix
  
  fdr_matrix <- data$fdr_value_matrix
  
  # Helper function to estimate space for the longest name, defined inside the main function
  estimate_space_for_text <- function(text_vector, rotation = 0, fontsize = 10) {
    # Find the length of the longest name
    longest_name <- max(nchar(text_vector))
    # Estimate space needed based on the length of the longest name
    width_in_inches <- longest_name * strwidth("W", units = "inches", gp = gpar(fontsize = fontsize))
    # Convert inches to cm and add a little extra space
    width_in_cm <- width_in_inches * 2.54 + 1
    return(unit(width_in_cm, "cm"))
  }
  
  # Estimate space needed for the longest row and column names
  row_names_margin <- estimate_space_for_text(rownames(chi2_matrix))
  column_names_margin <- estimate_space_for_text(colnames(chi2_matrix), rotation = 45)
  
  # Create color ramp
  pastel_colors <- colorRampPalette(c(palette))(100)
  
  # Define the heatmap with cell_fun to add text inside each cell
  ht <- Heatmap(chi2_matrix, name = "example", 
                col = pastel_colors,
                cluster_rows = T, cluster_columns = T,
                row_names_max_width = row_names_margin,
                column_names_max_height = column_names_margin,
                # heatmap_legend_param = list(direction = "horizontal"),
                heatmap_legend_param = list(
                  title = "log2(chi2 + 1)", # Add a title if needed
                  direction = "horizontal",
                  # title_position = "leftcenter", # Position the title to the left
                  legend_width = unit(6, "cm"), # Adjust the width as needed
                  legend_height = unit(1, "cm"), # Adjust the height if needed
                  labels_gp = gpar(fontsize = 10), # Adjust label font size if needed
                  title_gp = gpar(fontsize = 12) # Adjust title font size if needed
                ),
                # cell_fun = function(j, i, x, y, width, height, fill) {
                #   grid.text(sprintf("%.1f", chi2_matrix[i, j]), x, y, gp = gpar(fontsize = 10))
                # }
                cell_fun = function(j, i, x, y, width, height, fill) {
                  # Retrieve the value from the chi2_matrix and the number_overlap_matrix
                  chi2_value <- sprintf("%.1f", chi2_matrix[i, j])
                  overlap_value <- number_overlap_matrix[i, j]
                  
                  # Combine the chi2 value and the overlap value in the desired format
                  combined_text <- sprintf("%s (%d)", chi2_value, overlap_value)
                  
                  # Draw the combined text in the cell
                  grid.text(combined_text, x, y, gp = gpar(fontsize = 10))
                }
  )
  
  # Define your custom legends
  significance_legend <- Legend(
    labels = c("FDR < 0.05"),
    border = color_rect,
    direction = "horizontal",
    # Change the font size of the labels
    grid_width = unit(1, "cm"),
    # Change the width of the legend key
    grid_height = unit(1, "cm")
  )  
  
  # Now draw the heatmap with both legends
  draw(ht, 
       annotation_legend_side = "bottom",
       heatmap_legend_side = "bottom",
       merge_legend = TRUE,
       annotation_legend_list = list(significance_legend))
  
  # Draw the heatmap with the legend on the right
  # draw(ht, heatmap_legend_side = "left", annotation_legend_side = "bot")
  
  #############################################################################
  # Extract the order of rows and columns from the heatmap object
  row_order <- row_order(ht)
  col_order <- column_order(ht)
  
  # Use the order to reorder the original matrix
  fdr_matrix <- fdr_matrix[row_order, col_order]
  
  # Specific paired row and column names to highlight
  highlight_pairs <- 
    fdr_matrix %>% 
    as.data.frame() %>% 
    mutate(row_names = rownames(.)) %>% 
    gather(key = "col_names", value = "value", -row_names) %>% 
    filter(value < fdr_threshold) %>% 
    select(c(row_names, col_names))
  
  
  # Check if the number of rows and columns to be paired is the same
  if (nrow(highlight_pairs) != length(highlight_pairs$col_names)) {
    stop("The number of row names and column names to highlight must be the same.")
  }
  
  # Convert row and column names to indices
  row_indices <- match(highlight_pairs$row_names, rownames(fdr_matrix))
  col_indices <- match(highlight_pairs$col_names, colnames(fdr_matrix))
  
  # Filter out any NAs from unmatched names
  valid_pairs <- highlight_pairs[!is.na(row_indices) & !is.na(col_indices), ]
  valid_row_indices <- row_indices[!is.na(row_indices)]
  valid_col_indices <- col_indices[!is.na(col_indices)]
  
  ### code to indicate signifficant results
  # Loop over the pairs to draw the rectangles
  for (i in seq_along(valid_row_indices)) {
    row_index <- valid_row_indices[i]
    col_index <- valid_col_indices[i]
    
    # Calculate normalized positions
    x_center <- (col_index - 0.5) / ncol(fdr_matrix)
    y_center <- (nrow(fdr_matrix) - row_index + 0.5) / nrow(fdr_matrix)
    
    # Draw the rectangle
    decorate_heatmap_body("example", {
      grid.rect(
        x = x_center,
        y = y_center,
        width = 1 / ncol(fdr_matrix),
        height = 1 / nrow(fdr_matrix),
        just = "center",
        gp = gpar(col = color_rect, lwd = lwd_rect, alpha = alpha_rect, fill = NA)
      )
    })
  }
  
  ### code to filling of significant results
  # Draw the rectangle with the pattern fill if 'apply_filling' is TRUE
  if(apply_filling){
    # Loop over the pairs to draw the rectangles with patterns
    for (i in seq_along(valid_row_indices)) {
      row_index <- valid_row_indices[i]
      col_index <- valid_col_indices[i]
      
      # Calculate the coordinates
      x_left <- (col_index - 1) / ncol(fdr_matrix)
      y_bottom <- (nrow(fdr_matrix) - row_index) / nrow(fdr_matrix)
      width <- 1 / ncol(fdr_matrix)
      height <- 1 / nrow(fdr_matrix)
      
      # Draw the rectangle with the pattern fill
      decorate_heatmap_body("example", {
        # First draw a blank rectangle
        grid.rect(
          x = x_left + width / 2,
          y = y_bottom + height / 2,
          width = width,
          height = height,
          just = "center",
          gp = gpar(col = NA, fill = NA)
        )
        
        # Now draw small dots to create a dotted pattern
        dot_spacing <- width / 10 # Adjust the spacing of the dots here
        seq_x <- seq(x_left + dot_spacing/2, x_left + width, by = dot_spacing)
        seq_y <- seq(y_bottom + dot_spacing/2, y_bottom + height, by = dot_spacing)
        for (j in seq_along(seq_x)) {
          for (k in seq_along(seq_y)) {
            grid.points(
              x = seq_x[j],
              y = seq_y[k],
              pch = pch_filling, # Type of point, 16 is a filled circle
              size = unit(size_filling, "mm"), # Adjust the size of the dots here
              gp = gpar(col = color_filling, alpha = alpha_filling)
            )
          }
        }
      })
    }
  }
}