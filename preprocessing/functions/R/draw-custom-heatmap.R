draw_custom_heatmap <- function(data_list, 
                                data_type,
                                fdr_threshold = 0.05,
                                color_rect = "green",
                                fdr_thresholds = c(0.05, 0.01), # Now a vector of thresholds
                                color_rects = c("green", "red"), # Corresponding colors for each threshold
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
                                gene_list_sizes = T,
                                col_mapping_vector=NULL, 
                                row_mapping_vector=NULL,
                                col_significant=F,
                                overlap_threshold = 3,
                                ...
                                ) {
  
  if(gene_list_sizes){
    if(is.null(col_mapping_vector)) {
      
      data_list[[data_type]]$list$p_value_matrix %>% colnames -> col_mapping_vector
      names(col_mapping_vector) <- col_mapping_vector
      
      intersect(names(col_mapping_vector), names(data_list$gene_list_sizes)) %>%
        as.data.frame() %>%
        set_colnames("original_name") %>%
        mutate(name = col_mapping_vector[original_name],
               size = data_list$gene_list_sizes[original_name]) %>%
        mutate(new_name = paste0(name, " ", "(", size, ")   ")) %>%
        select(c(original_name, new_name)) -> col_mapping_df
      
      print(col_mapping_df)

      col_mapping_vector <- col_mapping_df$new_name
      names(col_mapping_vector) <- col_mapping_df$original_name
      
    } else {
      
      intersect(names(col_mapping_vector), names(data_list$gene_list_sizes)) %>%
        as.data.frame() %>%
        set_colnames("original_name") %>% 
        mutate(name = col_mapping_vector[original_name],
               size = data_list$gene_list_sizes[original_name]) %>% 
        mutate(new_name = paste0(name, " ", "(", size, ")   ")) %>% 
        select(c(original_name, new_name)) -> col_mapping_df
      
      col_mapping_vector <- col_mapping_df$new_name
      names(col_mapping_vector) <- col_mapping_df$original_name
      
    }
    
    if(is.null(row_mapping_vector)) {
      print("echo")
      
      data_list[[data_type]]$list$p_value_matrix %>% rownames() -> row_mapping_vector
      names(row_mapping_vector) <- row_mapping_vector
      
      print(row_mapping_vector)
      
      intersect(names(row_mapping_vector), names(data_list$gene_list_sizes)) %>%
        as.data.frame() %>%
        set_colnames("original_name") %>%
        mutate(name = row_mapping_vector[original_name],
               size = data_list$gene_list_sizes[original_name]) %>%
        mutate(new_name = paste0(name, " ", "(", size, ")")) %>%
        select(c(original_name, new_name)) -> row_mapping_df
      
      row_mapping_vector <- row_mapping_df$new_name
      names(row_mapping_vector) <- row_mapping_df$original_name
      
    } else {
      intersect(names(row_mapping_vector), names(data_list$gene_list_sizes)) %>%
        as.data.frame() %>%
        set_colnames("original_name") %>% 
        mutate(name = row_mapping_vector[original_name],
               size = data_list$gene_list_sizes[original_name]) %>% 
        # mutate(new_name = paste0(name, " ", "(", size, ")")) %>% 
        mutate(new_name = paste0(name, "    ", size, "    ")) %>% 
        select(c(original_name, new_name)) %>% 
        mutate(new_name = sapply(new_name, function(s) {
          # parts <- strsplit(s, "PMID", fixed = TRUE)[[1]]
          # result <- paste0("PMID", parts[2], ", ", parts[1])
          # Remove trailing comma
          result <- gsub(",\\s*$", "", s)
          # Add four spaces at the beginning
          result <- paste0("    ", result)
          result
        }))-> row_mapping_df
      
      row_mapping_vector <- row_mapping_df$new_name
      names(row_mapping_vector) <- row_mapping_df$original_name
      
    }
    
  }

  
  
  
    # Validate that the length of fdr_thresholds and color_rects are the same
  if (length(fdr_thresholds) != length(color_rects)) {
    stop("The number of thresholds and colors must be the same.")
  }
  
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
  pastel_colors <- colorRamp2(c(0, 2, 4, 6, 8), palette)
  
  
  # Define the heatmap with cell_fun to add text inside each cell
  ht <- Heatmap(chi2_matrix, name = "example", 
                col = pastel_colors,
                # cluster_rows = T, cluster_columns = T, 
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
                },
                ...
                )

  # Define custom legends for each threshold
  significance_legends <- lapply(seq_along(fdr_thresholds), function(t) {
    threshold <- fdr_thresholds[t]
    color <- color_rects[t]
    label <- paste0("FDR < ", threshold)
    
    Legend(
      labels = label,
      border = color,
      direction = "horizontal",
      # Change the font size of the labels
      grid_width = unit(1, "cm"),
      # Change the width of the legend key
      grid_height = unit(1, "cm")
    )
  })
  
  # Combine all legends into a single list
  all_legends <- do.call(c, significance_legends)
  
  # Now draw the heatmap with all legends
  draw(ht, 
       annotation_legend_side = "bottom",
       heatmap_legend_side = "bottom",
       merge_legend = TRUE,
       annotation_legend_list = all_legends)
  
  #############################################################################
  # Extract the order of rows and columns from the heatmap object
  row_order <- row_order(ht)
  col_order <- column_order(ht)
  
  # Use the order to reorder the original matrix
  fdr_matrix <- fdr_matrix[row_order, col_order]
  number_overlap_matrix <- number_overlap_matrix[row_order, col_order]
  
  # Loop over each threshold and its corresponding color
  for (t in seq_along(fdr_thresholds)) {
    fdr_threshold <- fdr_thresholds[t]
    color_rect <- color_rects[t]
    
    # Extract pairs that are significant at the current threshold
    highlight_pairs_fdr <- 
      fdr_matrix %>% 
      as.data.frame() %>% 
      mutate(row_names = rownames(.)) %>% 
      gather(key = "col_names", value = "value", -row_names) %>% 
      filter(value < fdr_threshold) %>% 
      select(c(row_names, col_names))
    
    highlight_pairs_number <- 
      number_overlap_matrix %>% 
      as.data.frame() %>% 
      mutate(row_names = rownames(.)) %>% 
      gather(key = "col_names", value = "value", -row_names) %>% 
      filter(value >= overlap_threshold) %>% 
      select(c(row_names, col_names))
    
    hihlight_pairs_inner <- inner_join(highlight_pairs_number, highlight_pairs_fdr)
    highlight_pairs <- hihlight_pairs_inner
    
    
    # Convert row and column names to indices
    row_indices <- match(highlight_pairs$row_names, rownames(fdr_matrix))
    col_indices <- match(highlight_pairs$col_names, colnames(fdr_matrix))
    
    # Filter out any NAs from unmatched names
    valid_pairs <- highlight_pairs[!is.na(row_indices) & !is.na(col_indices), ]
    valid_row_indices <- row_indices[!is.na(row_indices)]
    valid_col_indices <- col_indices[!is.na(col_indices)]
    
    
    # Loop over the pairs to draw the rectangles with the current color
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
  }
  

  
  ### code to filling of significant results
  # Draw the rectangle with the pattern fill if 'apply_filling' is TRUE
  fdr_threshold <- max(fdr_thresholds)
  
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
  
  decorate_heatmap_body("example", {
    grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 2))
    
  })
}

# 
# # # prepare data
# filter_matrices(chi2_results_phenotypes,
#                 !rownames(chi2_results_phenotypes$p_value_matrix) %in% tissues_clusters,
#                 tissues_clusters[10:27]) -> chi2_phenotypes_filter


# draw_custom_heatmap(
#   tmp,
#   data_type = "original_data",
#   col_mapping_vector = cluster_vector,
#   # row_mapping_vector =  setNames(papers_data_preprocessing$label2, papers_data_preprocessing$label),
#   fdr_threshold = 0.1,
#   fdr_thresholds = c(0.05, 0.0001),
#   color_rects =  c("green", "#FF00FF"),
#   color_rect = "green",
#   lwd_rect = 3,
#   alpha_rect = 1,
#   apply_filling = F,
#   color_filling = "gray",
#   alpha_filling = 0.6,
#   size_filling = 1,
#   pch_filling = 16,
#   col_significant = F,
#   gene_list_sizes = T
# )
draw_custom_heatmap_v2 <- function(data_list, 
                                   data_type,
                                   p_thresholds = c(0.05, 0.01),
                                   color_rects = c("green", "red"),
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
                                     "white"             = "#ffffcc",
                                     "pastel_orange"     = "#fed976",
                                     "pastel_red"        = "#fc4e2a"
                                   ),
                                   color_scale_range = NULL,
                                   gene_list_sizes = TRUE,
                                   col_mapping_vector = NULL, 
                                   row_mapping_vector = NULL,
                                   col_significant = FALSE,
                                   overlap_threshold = 3,
                                   ...
) {
  if (gene_list_sizes) {
    # kolumny
    if (is.null(col_mapping_vector)) {
      data_list[[data_type]]$list$p_value_matrix %>% colnames() -> col_mapping_vector
      names(col_mapping_vector) <- col_mapping_vector
      intersect(names(col_mapping_vector), names(data_list$gene_list_sizes)) %>%
        as.data.frame() %>%
        setNames("original_name") %>%
        mutate(name = col_mapping_vector[original_name],
               size = data_list$gene_list_sizes[original_name]) %>%
        mutate(new_name = paste0(name, " (", size, ")")) %>%
        select(original_name, new_name) -> col_mapping_df
      col_mapping_vector <- setNames(col_mapping_df$new_name, col_mapping_df$original_name)
    } else {
      intersect(names(col_mapping_vector), names(data_list$gene_list_sizes)) %>%
        as.data.frame() %>%
        setNames("original_name") %>%
        mutate(name = col_mapping_vector[original_name],
               size = data_list$gene_list_sizes[original_name]) %>%
        mutate(new_name = paste0(name, " (", size, ")")) %>%
        select(original_name, new_name) -> col_mapping_df
      col_mapping_vector <- setNames(col_mapping_df$new_name, col_mapping_df$original_name)
    }
    
    # wiersze
    if (is.null(row_mapping_vector)) {
      data_list[[data_type]]$list$p_value_matrix %>% rownames() -> row_mapping_vector
      names(row_mapping_vector) <- row_mapping_vector
      intersect(names(row_mapping_vector), names(data_list$gene_list_sizes)) %>%
        as.data.frame() %>%
        setNames("original_name") %>%
        mutate(name = row_mapping_vector[original_name],
               size = data_list$gene_list_sizes[original_name]) %>%
        mutate(new_name = paste0(name, " (", size, ")")) %>%
        select(original_name, new_name) -> row_mapping_df
      row_mapping_vector <- setNames(row_mapping_df$new_name, row_mapping_df$original_name)
    } else {
      intersect(names(row_mapping_vector), names(data_list$gene_list_sizes)) %>%
        as.data.frame() %>%
        setNames("original_name") %>%
        mutate(name = row_mapping_vector[original_name],
               size = data_list$gene_list_sizes[original_name]) %>%
        mutate(new_name = paste0(name, " (", size, ")")) %>%
        select(original_name, new_name) -> row_mapping_df
      row_mapping_vector <- setNames(row_mapping_df$new_name, row_mapping_df$original_name)
    }
  }
  
  if (length(p_thresholds) != length(color_rects)) {
    stop("The number of thresholds and colors must be the same.")
  }
  
  data <- data_list[[data_type]]$list
  cols_to_filter <- data_list[[data_type]]$cols
  rows_to_filter <- data_list[[data_type]]$rows
  
  if (col_significant) {
    data <- lapply(data, function(mat) mat[rows_to_filter, cols_to_filter])
  }
  
  replace_names_in_data(data, col_mapping_vector, row_mapping_vector) -> data
  
  chi2_matrix <- log2(data$chi2_value_matrix + 1)
  number_overlap_matrix <- data$number_overlap_matrix
  p_matrix <- data$p_value_matrix
  
  # kolorowanie - ustalanie zakresu skali
  if (is.null(color_scale_range)) {
    color_scale_range <- range(chi2_matrix, na.rm = TRUE)
  }
  
  n_colors <- length(palette)
  scale_values <- seq(color_scale_range[1], color_scale_range[2], length.out = n_colors)
  pastel_colors <- colorRamp2(scale_values, palette)
  
  row_names_margin <- unit(max(nchar(rownames(chi2_matrix))) * 0.25, "cm")
  column_names_margin <- unit(max(nchar(colnames(chi2_matrix))) * 0.25, "cm")
  
  ht <- Heatmap(
    chi2_matrix,
    name = "example",
    col = pastel_colors,
    row_names_max_width = row_names_margin,
    column_names_max_height = column_names_margin,
    heatmap_legend_param = list(
      title = "log2(chi2 + 1)",
      direction = "horizontal",
      legend_width = unit(6, "cm"),
      legend_height = unit(1, "cm"),
      labels_gp = gpar(fontsize = 10),
      title_gp = gpar(fontsize = 12)
    ),
    cell_fun = function(j, i, x, y, width, height, fill) {
      val <- sprintf("%.1f", chi2_matrix[i, j])
      ovl <- number_overlap_matrix[i, j]
      grid.text(sprintf("%s (%d)", val, ovl), x, y, gp = gpar(fontsize = 10))
    },
    ...
  )
  
  legends <- lapply(seq_along(p_thresholds), function(t) {
    Legend(
      labels = paste0("p < ", p_thresholds[t]),
      border = color_rects[t],
      direction = "horizontal",
      grid_width = unit(1, "cm"),
      grid_height = unit(1, "cm")
    )
  })
  
  draw(ht,
       annotation_legend_side = "bottom",
       heatmap_legend_side = "bottom",
       merge_legend = TRUE,
       annotation_legend_list = do.call(c, legends))
  
  row_order <- row_order(ht)
  col_order <- column_order(ht)
  p_matrix <- p_matrix[row_order, col_order]
  number_overlap_matrix <- number_overlap_matrix[row_order, col_order]
  
  for (t in seq_along(p_thresholds)) {
    threshold <- p_thresholds[t]
    color <- color_rects[t]
    
    highlight_pairs_p <- as.data.frame(p_matrix) %>%
      tibble::rownames_to_column("row_names") %>%
      tidyr::pivot_longer(-row_names, names_to = "col_names", values_to = "value") %>%
      filter(value < threshold)
    
    highlight_pairs_overlap <- as.data.frame(number_overlap_matrix) %>%
      tibble::rownames_to_column("row_names") %>%
      tidyr::pivot_longer(-row_names, names_to = "col_names", values_to = "value") %>%
      filter(value >= overlap_threshold)
    
    highlight_pairs <- inner_join(highlight_pairs_p, highlight_pairs_overlap, by = c("row_names", "col_names"))
    
    row_idx <- match(highlight_pairs$row_names, rownames(p_matrix))
    col_idx <- match(highlight_pairs$col_names, colnames(p_matrix))
    
    for (i in seq_along(row_idx)) {
      x <- (col_idx[i] - 0.5) / ncol(p_matrix)
      y <- (nrow(p_matrix) - row_idx[i] + 0.5) / nrow(p_matrix)
      
      decorate_heatmap_body("example", {
        grid.rect(x = x, y = y,
                  width = 1 / ncol(p_matrix),
                  height = 1 / nrow(p_matrix),
                  just = "center",
                  gp = gpar(col = color, lwd = lwd_rect, alpha = alpha_rect, fill = NA))
      })
    }
  }
  
  if (apply_filling && exists("highlight_pairs")) {
    row_idx <- match(highlight_pairs$row_names, rownames(p_matrix))
    col_idx <- match(highlight_pairs$col_names, colnames(p_matrix))
    
    for (i in seq_along(row_idx)) {
      x_left <- (col_idx[i] - 1) / ncol(p_matrix)
      y_bottom <- (nrow(p_matrix) - row_idx[i]) / nrow(p_matrix)
      width <- 1 / ncol(p_matrix)
      height <- 1 / nrow(p_matrix)
      
      decorate_heatmap_body("example", {
        grid.rect(x = x_left + width / 2, y = y_bottom + height / 2,
                  width = width, height = height, just = "center",
                  gp = gpar(col = NA, fill = NA))
        
        seq_x <- seq(x_left + width / 20, x_left + width, by = width / 10)
        seq_y <- seq(y_bottom + height / 20, y_bottom + height, by = height / 10)
        for (xx in seq_x) {
          for (yy in seq_y) {
            grid.points(x = xx, y = yy,
                        pch = pch_filling,
                        size = unit(size_filling, "mm"),
                        gp = gpar(col = color_filling, alpha = alpha_filling))
          }
        }
      })
    }
  }
  
  decorate_heatmap_body("example", {
    grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 2))
  })
}

draw_custom_heatmap_v3 <- function(data_list, 
                                   data_type,
                                   p_thresholds = c(0.05, 0.01),
                                   color_rects = c("green", "red"),
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
                                     "white"             = "#ffffcc",
                                     "pastel_orange"     = "#fed976",
                                     "pastel_red"        = "#fc4e2a"
                                   ),
                                   color_scale_range = NULL,
                                   gene_list_sizes = TRUE,
                                   col_only_n_genes = FALSE,
                                   col_mapping_vector = NULL, 
                                   row_mapping_vector = NULL,
                                   col_significant = FALSE,
                                   overlap_threshold = 3,
                                   ...) {
  if (gene_list_sizes) {
    # kolumny
    data_cols <- colnames(data_list[[data_type]]$list$p_value_matrix)
    if (is.null(col_mapping_vector)) {
      col_mapping_vector <- setNames(data_cols, data_cols)
    }
    common_cols <- intersect(names(col_mapping_vector), names(data_list$gene_list_sizes))
    col_mapping_df <- data.frame(original_name = common_cols) %>%
      mutate(
        name = col_mapping_vector[original_name],
        size = data_list$gene_list_sizes[original_name],
        new_name = if (col_only_n_genes) {
          as.character(size)
        } else {
          paste0(name, " (", size, ")")
        }
      )
    col_mapping_vector <- setNames(col_mapping_df$new_name, col_mapping_df$original_name)
    
    # wiersze
    data_rows <- rownames(data_list[[data_type]]$list$p_value_matrix)
    if (is.null(row_mapping_vector)) {
      row_mapping_vector <- setNames(data_rows, data_rows)
    }
    common_rows <- intersect(names(row_mapping_vector), names(data_list$gene_list_sizes))
    row_mapping_df <- data.frame(original_name = common_rows) %>%
      mutate(
        name = row_mapping_vector[original_name],
        size = data_list$gene_list_sizes[original_name],
        new_name = paste0(name, " (", size, ")")
      )
    row_mapping_vector <- setNames(row_mapping_df$new_name, row_mapping_df$original_name)
  }
  
  if (length(p_thresholds) != length(color_rects)) {
    stop("The number of thresholds and colors must be the same.")
  }
  
  data <- data_list[[data_type]]$list
  cols_to_filter <- data_list[[data_type]]$cols
  rows_to_filter <- data_list[[data_type]]$rows
  
  if (col_significant) {
    data <- lapply(data, function(mat) mat[rows_to_filter, cols_to_filter])
  }
  
  replace_names_in_data(data, col_mapping_vector, row_mapping_vector) -> data
  
  chi2_matrix <- log2(data$chi2_value_matrix + 1)
  number_overlap_matrix <- data$number_overlap_matrix
  p_matrix <- data$p_value_matrix
  
  if (is.null(color_scale_range)) {
    color_scale_range <- range(chi2_matrix, na.rm = TRUE)
  }
  
  n_colors <- length(palette)
  scale_values <- seq(color_scale_range[1], color_scale_range[2], length.out = n_colors)
  pastel_colors <- colorRamp2(scale_values, palette)
  
  row_names_margin <- unit(max(nchar(rownames(chi2_matrix))) * 0.25, "cm")
  column_names_margin <- unit(max(nchar(colnames(chi2_matrix))) * 0.25, "cm")
  
  ht <- Heatmap(
    chi2_matrix,
    name = "example",
    col = pastel_colors,
    row_names_max_width = row_names_margin,
    column_names_max_height = column_names_margin,
    heatmap_legend_param = list(
      title = "log2(chi2 + 1)",
      direction = "horizontal",
      legend_width = unit(6, "cm"),
      legend_height = unit(1, "cm"),
      labels_gp = gpar(fontsize = 10),
      title_gp = gpar(fontsize = 12)
    ),
    cell_fun = function(j, i, x, y, width, height, fill) {
      val <- sprintf("%.1f", chi2_matrix[i, j])
      ovl <- number_overlap_matrix[i, j]
      grid.text(sprintf("%s (%d)", val, ovl), x, y, gp = gpar(fontsize = 10))
    },
    ...
  )
  
  legends <- lapply(seq_along(p_thresholds), function(t) {
    Legend(
      labels = paste0("p < ", p_thresholds[t]),
      border = color_rects[t],
      direction = "horizontal",
      grid_width = unit(1, "cm"),
      grid_height = unit(1, "cm")
    )
  })
  
  draw(ht,
       annotation_legend_side = "bottom",
       heatmap_legend_side = "bottom",
       merge_legend = TRUE,
       annotation_legend_list = do.call(c, legends))
  
  row_order <- row_order(ht)
  col_order <- column_order(ht)
  p_matrix <- p_matrix[row_order, col_order]
  number_overlap_matrix <- number_overlap_matrix[row_order, col_order]
  
  for (t in seq_along(p_thresholds)) {
    threshold <- p_thresholds[t]
    color <- color_rects[t]
    
    highlight_pairs_p <- as.data.frame(p_matrix) %>%
      tibble::rownames_to_column("row_names") %>%
      tidyr::pivot_longer(-row_names, names_to = "col_names", values_to = "value") %>%
      filter(value < threshold)
    
    highlight_pairs_overlap <- as.data.frame(number_overlap_matrix) %>%
      tibble::rownames_to_column("row_names") %>%
      tidyr::pivot_longer(-row_names, names_to = "col_names", values_to = "value") %>%
      filter(value >= overlap_threshold)
    
    highlight_pairs <- inner_join(highlight_pairs_p, highlight_pairs_overlap, by = c("row_names", "col_names"))
    
    row_idx <- match(highlight_pairs$row_names, rownames(p_matrix))
    col_idx <- match(highlight_pairs$col_names, colnames(p_matrix))
    
    for (i in seq_along(row_idx)) {
      x <- (col_idx[i] - 0.5) / ncol(p_matrix)
      y <- (nrow(p_matrix) - row_idx[i] + 0.5) / nrow(p_matrix)
      
      decorate_heatmap_body("example", {
        grid.rect(x = x, y = y,
                  width = 1 / ncol(p_matrix),
                  height = 1 / nrow(p_matrix),
                  just = "center",
                  gp = gpar(col = color, lwd = lwd_rect, alpha = alpha_rect, fill = NA))
      })
    }
  }
  
  if (apply_filling && exists("highlight_pairs")) {
    row_idx <- match(highlight_pairs$row_names, rownames(p_matrix))
    col_idx <- match(highlight_pairs$col_names, colnames(p_matrix))
    
    for (i in seq_along(row_idx)) {
      x_left <- (col_idx[i] - 1) / ncol(p_matrix)
      y_bottom <- (nrow(p_matrix) - row_idx[i]) / nrow(p_matrix)
      width <- 1 / ncol(p_matrix)
      height <- 1 / nrow(p_matrix)
      
      decorate_heatmap_body("example", {
        grid.rect(x = x_left + width / 2, y = y_bottom + height / 2,
                  width = width, height = height, just = "center",
                  gp = gpar(col = NA, fill = NA))
        
        seq_x <- seq(x_left + width / 20, x_left + width, by = width / 10)
        seq_y <- seq(y_bottom + height / 20, y_bottom + height, by = height / 10)
        for (xx in seq_x) {
          for (yy in seq_y) {
            grid.points(x = xx, y = yy,
                        pch = pch_filling,
                        size = unit(size_filling, "mm"),
                        gp = gpar(col = color_filling, alpha = alpha_filling))
          }
        }
      })
    }
  }
  
  decorate_heatmap_body("example", {
    grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 2))
  })
}
