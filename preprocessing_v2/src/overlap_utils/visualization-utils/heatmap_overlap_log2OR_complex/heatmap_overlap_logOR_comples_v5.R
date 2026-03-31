heatmap_overlap_log2OR_complex <- function(
    data_list,
    data_type,
    
    rows_to_filter = NULL,
    cols_to_filter = NULL,
    
    row_mapper = NULL,
    col_mapper = NULL,
    
    row_order_original = NULL,
    col_order_original = NULL,
    
    # ---- choose matrix for COLORS (and LABEL VALUES) by key in data_list[[data_type]]$list ----
    color_key = "log2_odds_ratio_matrix",   # e.g. "log_chi2_value_matrix", "p_value_matrix", "number_overlap_matrix"
    custom_color_matrix = NULL,
    color_legend_title = NULL,
    color_symmetric = TRUE,
    
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
    
    row_dend_height = grid::unit(20, "mm"),
    col_dend_height = grid::unit(20, "mm"),
    row_names_width = grid::unit(80, "mm"),
    col_names_height = grid::unit(50, "mm"),
    
    tile_gap = 1,
    corner_size = 0.28,
    
    save_to_svg = NULL,
    force_create_directory = FALSE,
    svg_width = 18,
    svg_height = 14
){
  
  triangle_mode <- match.arg(triangle_mode)
  
  # ---- SAVE ----
  if (!is.null(save_to_svg)) {
    dir_path <- dirname(save_to_svg)
    if (!dir.exists(dir_path)) {
      if (force_create_directory) dir.create(dir_path, recursive = TRUE)
      else stop("Directory does not exist: ", dir_path,
                "\nSet force_create_directory = TRUE to create it.")
    }
    message("📁 Saving output to SVG: ", save_to_svg)
    grDevices::svg(save_to_svg, width = svg_width, height = svg_height)
  }
  
  # ---- LOAD DATA ----
  data <- data_list[[data_type]]$list
  
  log2or     <- as.matrix(data$log2_odds_ratio_matrix)
  overlap    <- as.matrix(data$number_overlap_matrix)
  pval       <- as.matrix(data$p_value_matrix)
  gene_sizes <- data_list$gene_list_sizes
  
  # ---- GET COLOR MATRIX EARLY (BEFORE FILTERS/MAPPING) ----
  if (identical(color_key, "custom")) {
    if (is.null(custom_color_matrix)) stop("color_key='custom' but custom_color_matrix is NULL.")
    hm_mat <- as.matrix(custom_color_matrix)
  } else {
    if (!color_key %in% names(data)) {
      stop("color_key not found in data_list[[data_type]]$list: ", color_key,
           "\nAvailable keys:\n- ", paste(names(data), collapse = "\n- "))
    }
    hm_mat <- as.matrix(data[[color_key]])
  }
  
  # ---- BASIC DIM CHECK (raw, before filtering) ----
  if (!all(dim(hm_mat) == dim(log2or))) {
    stop("Selected color matrix has different dimensions than log2_odds_ratio_matrix.\n",
         "color_key=", color_key, "\n",
         "dim(color)=", paste(dim(hm_mat), collapse="x"),
         " vs dim(log2or)=", paste(dim(log2or), collapse="x"))
  }
  
  # ---- FILTER ROWS / COLS (apply to ALL matrices incl hm_mat) ----
  if (!is.null(rows_to_filter)) {
    keep <- intersect(rownames(log2or), rows_to_filter)
    log2or  <- log2or[keep,,drop=FALSE]
    overlap <- overlap[keep,,drop=FALSE]
    pval    <- pval[keep,,drop=FALSE]
    hm_mat  <- hm_mat[keep,,drop=FALSE]
  }
  
  if (!is.null(cols_to_filter)) {
    keep <- intersect(colnames(log2or), cols_to_filter)
    log2or  <- log2or[,keep,drop=FALSE]
    overlap <- overlap[,keep,drop=FALSE]
    pval    <- pval[,keep,drop=FALSE]
    hm_mat  <- hm_mat[,keep,drop=FALSE]
  }
  
  # ---- NAME MAPPING ----
  current_rows <- rownames(log2or)
  current_cols <- colnames(log2or)
  
  safe_lookup <- function(x) gene_sizes[match(x, names(gene_sizes))]
  row_sizes <- safe_lookup(current_rows)
  col_sizes <- safe_lookup(current_cols)
  
  mapped_rows <- if (!is.null(row_mapper)) {
    ifelse(current_rows %in% names(row_mapper), row_mapper[current_rows], current_rows)
  } else current_rows
  
  mapped_cols <- if (!is.null(col_mapper)) {
    ifelse(current_cols %in% names(col_mapper), col_mapper[current_cols], current_cols)
  } else current_cols
  
  row_labels <- paste0(mapped_rows, " (", row_sizes, ")")
  col_labels <- paste0(mapped_cols, " (", col_sizes, ")")
  
  # set the same labels everywhere
  rownames(log2or)  <- row_labels; colnames(log2or)  <- col_labels
  rownames(overlap) <- row_labels; colnames(overlap) <- col_labels
  rownames(pval)    <- row_labels; colnames(pval)    <- col_labels
  rownames(hm_mat)  <- row_labels; colnames(hm_mat)  <- col_labels
  
  # ---- MANUAL ORDER (now works because labels already updated) ----
  if (!is.null(row_order_original)) {
    original_to_full_rows <- setNames(row_labels, mapped_rows)
    desired_rows <- original_to_full_rows[row_order_original]
    if (any(is.na(desired_rows))) {
      stop("Not all names from row_order_original found after mapping:\n",
           paste(row_order_original[is.na(desired_rows)], collapse=", "))
    }
    log2or  <- log2or[desired_rows,,drop=FALSE]
    overlap <- overlap[desired_rows,,drop=FALSE]
    pval    <- pval[desired_rows,,drop=FALSE]
    hm_mat  <- hm_mat[desired_rows,,drop=FALSE]
  }
  
  if (!is.null(col_order_original)) {
    original_to_full_cols <- setNames(col_labels, mapped_cols)
    desired_cols <- original_to_full_cols[col_order_original]
    if (any(is.na(desired_cols))) {
      stop("Not all names from col_order_original found after mapping:\n",
           paste(col_order_original[is.na(desired_cols)], collapse=", "))
    }
    log2or  <- log2or[,desired_cols,drop=FALSE]
    overlap <- overlap[,desired_cols,drop=FALSE]
    pval    <- pval[,desired_cols,drop=FALSE]
    hm_mat  <- hm_mat[,desired_cols,drop=FALSE]
  }
  
  # ---- TRIANGLE MODE ----
  if (triangle_mode != "full" && nrow(log2or) == ncol(log2or)) {
    n <- nrow(log2or)
    mask <- outer(seq_len(n), seq_len(n), function(i,j){
      if (triangle_mode == "upper") i >= j else i <= j
    })
    log2or[mask]  <- NA
    overlap[mask] <- NA
    pval[mask]    <- NA
    hm_mat[mask]  <- NA
  }
  
  # ---- LABEL MATRIX: ALWAYS value from color_key + (overlap) ----
  if (identical(color_key, "number_overlap_matrix")) {
    label_mat <- ifelse(!is.na(overlap), sprintf("%d", overlap), "")
    
  } else if (identical(color_key, "p_value_matrix")) {
    # p-values w notacji naukowej
    label_mat <- ifelse(!is.na(hm_mat), sprintf("%.1e (%d)", hm_mat, overlap), "")
    
  } else {
    label_mat <- ifelse(!is.na(hm_mat), sprintf("%.2f (%d)", hm_mat, overlap), "")
  }
  
  # ---- TEXT COLOR MATRIX (contrast uses hm_mat) ----
  text_col_mat <- matrix(text_color, nrow=nrow(hm_mat), ncol=ncol(hm_mat))
  if (!is.null(text_contrast_range)) {
    inside <- hm_mat >= text_contrast_range[1] & hm_mat <= text_contrast_range[2]
    text_col_mat[inside] <- "black"
    text_col_mat[!inside & !is.na(hm_mat)] <- "white"
  }
  
  # ---- SIGNIFICANCE ----
  sig_list <- lapply(seq_along(p_thresholds), function(i){
    (pval < p_thresholds[i]) & (overlap >= overlap_threshold) & !is.na(pval)
  })
  
  # ---- COLOR SCALE (based on hm_mat) ----
  if (is.null(color_scale_range)) {
    if (isTRUE(color_symmetric)) {
      max_abs <- max(abs(hm_mat), na.rm=TRUE)
      color_scale_range <- c(-max_abs, max_abs)
    } else {
      color_scale_range <- range(hm_mat, na.rm=TRUE)
    }
  }
  breaks <- seq(color_scale_range[1], color_scale_range[2], length.out=length(palette))
  col_fun <- circlize::colorRamp2(breaks, palette)
  
  if (is.null(color_legend_title)) {
    color_legend_title <- if (identical(color_key, "log2_odds_ratio_matrix")) "log2(OR)" else color_key
  }
  
  # ---- HEATMAP ----
  ht <- ComplexHeatmap::Heatmap(
    hm_mat,
    name = color_legend_title,
    col = col_fun,
    na_col = NA,
    border = TRUE,
    border_gp = grid::gpar(col="black", lwd=2),
    
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
    row_names_gp = grid::gpar(fontsize=text_size_axis),
    column_names_side = "top",
    column_names_gp = grid::gpar(fontsize=text_size_axis, rot=axis_text_angle, just=axis_text_hjust),
    
    column_title = title,
    column_title_gp = grid::gpar(fontsize = title_size),
    
    heatmap_legend_param = list(
      title = color_legend_title,
      direction = "horizontal",
      legend_height = grid::unit(6, "mm"),
      title_gp = grid::gpar(fontsize=text_size_legend),
      labels_gp = grid::gpar(fontsize=text_size_legend)
    ),
    
    rect_gp = grid::gpar(col="white", lwd=tile_gap),
    
    cell_fun = function(j,i,x,y,w,h,fill){
      
      if (is.na(hm_mat[i,j])) return()
      
      grid::grid.text(
        label_mat[i,j], x, y,
        gp = grid::gpar(col=text_col_mat[i,j], fontsize=text_size_tile * 2.5)
      )
      
      for (k in seq_along(sig_list)) {
        if (sig_list[[k]][i,j]) {
          grid::grid.rect(
            x, y, width=w, height=h,
            gp=grid::gpar(col=color_rects[k], lwd=rect_lwd, fill=NA)
          )
          grid::grid.polygon(
            x = grid::unit.c(x - w/2, x - w/2, x - w/2 + w*corner_size),
            y = grid::unit.c(y + h/2, y + h/2 - h*corner_size, y + h/2),
            gp = grid::gpar(fill=color_rects[k], col=color_rects[k], lwd=1)
          )
        }
      }
    }
  )
  
  legend_rectangles <- ComplexHeatmap::Legend(
    at = paste0("p < ", p_thresholds),
    type = "points",
    pch = 22,
    size = grid::unit(6,"mm"),
    nrow = 1, ncol = 2,
    direction = "horizontal",
    legend_gp = grid::gpar(fill = NA, col = color_rects, lwd = rect_lwd),
    title = "Significance",
    title_gp = grid::gpar(fontsize = text_size_legend),
    labels_gp = grid::gpar(fontsize = text_size_legend),
    title_position = "leftcenter"
  )
  
  ComplexHeatmap::draw(
    ht,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom",
    annotation_legend_list = list(legend_rectangles),
    padding = grid::unit(c(10,10,30,40),"mm")
  )
  
  if (!is.null(save_to_svg)) {
    grDevices::dev.off()
    message("✔ SVG saved: ", save_to_svg)
  }
  
  invisible(ht)
}
