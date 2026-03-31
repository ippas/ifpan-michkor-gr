library(ComplexHeatmap)
library(circlize)
library(grid)

heatmap_overlap_log2OR_complex <- function(
    data_list,
    data_type,
    
    rows_to_filter = NULL,
    cols_to_filter = NULL,
    
    row_mapper = NULL,
    col_mapper = NULL,
    
    p_thresholds = c(0.05, 0.01),
    color_rects = c("#97C426", "#2F4603"),
    overlap_threshold = 3,
    rect_lwd = 2.5,
    
    triangle_mode = c("full","upper","lower"),
    
    palette = c("#07243e", "white", "#8b0000"),
    color_scale_range = NULL,
    
    text_contrast_range = NULL,
    text_color = "black",
    text_size_tile = 3,
    text_size_axis = 10,
    text_size_legend = 10,
    axis_text_angle = 45,
    axis_text_hjust = 0,
    
    title = NULL,
    title_size = 14,
    
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_dendrograms = TRUE,
    cluster_method = "complete",
    
    # dendrogram
    row_dend_height = unit(20, "mm"),
    col_dend_height = unit(20, "mm"),
    row_names_width = unit(80, "mm"),
    col_names_height = unit(50, "mm"),
    
    tile_gap = 1,
    
    # nowy parametr rogu
    corner_size = 0.28,
    
    save_to_svg = NULL,
    force_create_directory = FALSE,
    svg_width = 18,
    svg_height = 14
){
  
  triangle_mode <- match.arg(triangle_mode)
  
  # ZAPIS
  if (!is.null(save_to_svg)) {
    dir_path <- dirname(save_to_svg)
    if (!dir.exists(dir_path)) {
      if (force_create_directory) {
        dir.create(dir_path, recursive = TRUE)
      } else {
        stop("Directory does not exist: ", dir_path,
             "\nSet force_create_directory = TRUE to create it.")
      }
    }
    message("📁 Saving output to SVG: ", save_to_svg)
    svg(save_to_svg, width = svg_width, height = svg_height)
  }
  
  # DANE
  data <- data_list[[data_type]]$list
  log2or   <- as.matrix(data$log2_odds_ratio_matrix)
  overlap  <- as.matrix(data$number_overlap_matrix)
  pval     <- as.matrix(data$p_value_matrix)
  gene_sizes <- data_list$gene_list_sizes
  
  # FILTRY
  if (!is.null(rows_to_filter)) {
    keep <- intersect(rownames(log2or), rows_to_filter)
    log2or <- log2or[keep,,drop=FALSE]
    overlap <- overlap[keep,,drop=FALSE]
    pval <- pval[keep,,drop=FALSE]
  }
  
  if (!is.null(cols_to_filter)) {
    keep <- intersect(colnames(log2or), cols_to_filter)
    log2or <- log2or[,keep,drop=FALSE]
    overlap <- overlap[,keep,drop=FALSE]
    pval <- pval[,keep,drop=FALSE]
  }
  
  # MAPOWANIE NAZW
  current_rows <- rownames(log2or)
  current_cols <- colnames(log2or)
  
  safe_lookup <- function(x) gene_sizes[match(x, names(gene_sizes))]
  row_sizes <- safe_lookup(current_rows)
  col_sizes <- safe_lookup(current_cols)
  
  mapped_rows <- ifelse(current_rows %in% names(row_mapper),
                        row_mapper[current_rows],
                        current_rows)
  mapped_cols <- ifelse(current_cols %in% names(col_mapper),
                        col_mapper[current_cols],
                        current_cols)
  
  row_labels <- paste0(mapped_rows, " (", row_sizes, ")")
  col_labels <- paste0(mapped_cols, " (", col_sizes, ")")
  
  rownames(log2or) <- row_labels
  colnames(log2or) <- col_labels
  rownames(overlap) <- row_labels
  colnames(overlap) <- col_labels
  rownames(pval) <- row_labels
  colnames(pval) <- col_labels
  
  # TRIANGLE MODE
  if (triangle_mode != "full" && nrow(log2or) == ncol(log2or)) {
    n <- nrow(log2or)
    mask <- outer(seq_len(n), seq_len(n), function(i,j){
      if (triangle_mode == "upper") i >= j else i <= j
    })
    log2or[mask] <- NA
    overlap[mask] <- NA
    pval[mask] <- NA
  }
  
  # LABELS
  label_mat <- ifelse(!is.na(log2or),
                      sprintf("%.2f (%d)",log2or,overlap),"")
  
  text_col_mat <- matrix(text_color, nrow=nrow(log2or), ncol=ncol(log2or))
  if (!is.null(text_contrast_range)) {
    inside <- log2or >= text_contrast_range[1] &
      log2or <= text_contrast_range[2]
    text_col_mat[inside] <- "black"
    text_col_mat[!inside & !is.na(log2or)] <- "white"
  }
  
  # SIGNIFICANCE
  sig_list <- lapply(seq_along(p_thresholds), function(i){
    (pval < p_thresholds[i]) & (overlap >= overlap_threshold) & !is.na(pval)
  })
  
  # KOLORY
  if (is.null(color_scale_range)) {
    max_abs <- max(abs(log2or),na.rm=TRUE)
    color_scale_range <- c(-max_abs, max_abs)
  }
  
  breaks <- seq(color_scale_range[1], color_scale_range[2], length.out=length(palette))
  col_fun <- circlize::colorRamp2(breaks, palette)
  
  # HEATMAP
  ht <- Heatmap(
    log2or,
    name="log2(OR)",
    col=col_fun,
    na_col=NA,
    border = TRUE,
    border_gp = gpar(col="black", lwd=2),
    cluster_rows = cluster_rows,
    cluster_columns = cluster_cols,
    clustering_method_rows = cluster_method,
    clustering_method_columns = cluster_method,
    show_row_dend = show_dendrograms,
    show_column_dend = show_dendrograms,
    row_dend_width = row_dend_height,
    column_dend_height = col_dend_height,
    row_names_max_width = row_names_width,
    column_names_max_height = col_names_height,
    row_names_side = "right",
    row_names_gp = gpar(fontsize=text_size_axis),
    column_names_side = "top",
    column_names_gp = gpar(fontsize=text_size_axis,
                           rot=axis_text_angle,
                           just=axis_text_hjust),
    column_title = title,
    column_title_gp = gpar(fontsize = title_size),
    heatmap_legend_param = list(
      title = "log2(OR)",
      direction = "horizontal",
      legend_height = unit(6, "mm"),
      title_gp = gpar(fontsize=text_size_legend),
      labels_gp = gpar(fontsize=text_size_legend)
    ),
    rect_gp = gpar(col="white", lwd=tile_gap),
    
    # ------------------------
    #  CELL_DRAWER + RÓG
    # ------------------------
    cell_fun = function(j,i,x,y,w,h,fill){
      if (is.na(log2or[i,j])) return()
      
      grid.text(label_mat[i,j], x,y,
                gp=gpar(col=text_col_mat[i,j],
                        fontsize=text_size_tile * 2.5))
      
      for (k in seq_along(sig_list)) {
        if (sig_list[[k]][i,j]) {
          # standardowa ramka
          grid.rect(x,y,width=w,height=h,
                    gp=gpar(col=color_rects[k], lwd=rect_lwd, fill=NA))
          
          # 🔥 NOWY ELEMENT - DOG EAR CORNER
          grid.polygon(
            x = unit.c(x - w/2, x - w/2, x - w/2 + w*corner_size),
            y = unit.c(y + h/2, y + h/2 - h*corner_size, y + h/2),
            gp = gpar(fill=color_rects[k], col=color_rects[k], lwd=1)
          )
        }
      }
    }
  )
  
  legend_rectangles <- Legend(
    at = paste0("p < ", p_thresholds),
    type = "points",
    pch = 22,
    size = unit(6,"mm"),
    nrow = 1, ncol = 2,
    direction = "horizontal",
    legend_gp = gpar(fill = NA, col = color_rects, lwd = rect_lwd),
    title = "Significance",
    title_gp = gpar(fontsize = text_size_legend),
    labels_gp = gpar(fontsize = text_size_legend),
    title_position = "leftcenter"
  )
  
  draw(ht,
       heatmap_legend_side = "bottom",
       annotation_legend_side = "bottom",
       annotation_legend_list = list(legend_rectangles),
       padding = unit(c(10,10,30,40),"mm"))
  
  if (!is.null(save_to_svg)) {
    dev.off()
    message("✔ SVG saved: ", save_to_svg)
  }
}
