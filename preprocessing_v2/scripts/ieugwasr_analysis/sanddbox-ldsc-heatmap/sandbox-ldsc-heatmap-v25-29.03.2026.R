create_rg_heatmap <- function(
    selected_traits,
    ldsc_results_annotated,
    metadata_ieu_EURsampleSize10000,
    
    title = "Genetic correlation heatmap (LDSC)",
    
    label_size = 8.5,
    rg_size = 7.5,
    star_size = 9,
    
    show_rg_values = TRUE,
    show_significance_stars = TRUE,
    
    rg_text_y_offset_mm = 1.0,
    star_text_y_offset_mm = -1.8,
    
    low_color  = "#2166AC",
    mid_color  = "#F7F7F7",
    high_color = "#B2182B",
    absolute_high_color = "#B2182B",
    
    row_labels_side = "left",
    col_labels_side = "top",
    
    cluster = TRUE,
    distance_method = "euclidean",
    clustering_method = "complete",
    
    absolute_correlation = FALSE,
    clustering_input = c("signed_rg", "absolute_rg", "displayed_matrix"),
    
    legend_gap_mm = 10,
    
    top_padding_mm = 8,
    right_padding_mm = 8,
    bottom_padding_mm = 8,
    left_padding_mm = 8,
    
    output_svg = NULL,
    svg_width = NULL,
    svg_height = NULL,
    
    output_png = NULL,
    png_width = NULL,
    png_height = NULL,
    png_res = 600,
    
    text_white_high_threshold = NULL,
    text_white_low_threshold = NULL,
    
    row_dend_width_mm = 28,
    column_dend_height_mm = 24,
    
    color_row_dend_branches = FALSE,
    color_col_dend_branches = FALSE,
    
    row_dend_k = NULL,
    col_dend_k = NULL,
    
    row_dend_colors = NULL,
    col_dend_colors = NULL,
    
    use_row_cluster_gaps = FALSE,
    use_col_cluster_gaps = FALSE,
    
    row_gap_color = "black",
    col_gap_color = "black",
    
    row_gap_width_mm = 1.2,
    col_gap_width_mm = 1.2,
    
    row_names_max_width_mm_cap = 85,
    column_names_max_height_mm_cap = 70,
    
    label_fields_rows = c("id", "trait"),
    label_fields_cols = c("id", "trait"),
    label_separator = " | ",
    
    verbose = TRUE,
    return_data = FALSE
) {
  
  # =========================================
  # 0. Helper for messages
  # =========================================
  .total_steps <- 8L
  .current_step <- 0L
  .last_step_time <- Sys.time()
  
  log_step <- function(msg) {
    if (isTRUE(verbose)) {
      now <- Sys.time()
      elapsed <- round(as.numeric(difftime(now, .last_step_time, units = "secs")), 2)
      .current_step <<- .current_step + 1L
      
      cat(
        sprintf(
          "[%s] [%d/%d] %s (%.2fs)\n",
          format(now, "%H:%M:%S"),
          .current_step,
          .total_steps,
          msg,
          elapsed
        )
      )
      
      .last_step_time <<- now
    }
  }
  
  # =========================================
  # 1. Packages + validation
  # =========================================
  suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(circlize)
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(grid)
    library(dendextend)
  })
  
  available_distance_methods <- c(
    "euclidean", "maximum", "manhattan",
    "canberra", "binary", "minkowski"
  )
  
  available_clustering_methods <- c(
    "complete", "single", "average", "mcquitty",
    "median", "centroid", "ward.D", "ward.D2"
  )
  
  available_label_fields <- c(
    "id",
    "trait",
    "category",
    "subcategory",
    "nsnp",
    "sample_size",
    "coverage",
    "year",
    "author",
    "sex",
    "population",
    "unit",
    "ncase",
    "ncontrol",
    "consortium",
    "study_design",
    "pmid",
    "doi",
    "ontology",
    "group_name",
    "build",
    "covariates",
    "priority",
    "mr",
    "qc_prior_to_upload",
    "note",
    "sd"
  )
  
  valid_svg_units <- c("in", "cm", "mm", "px")
  valid_png_units <- c("in", "cm", "mm", "px")
  
  clustering_input <- match.arg(clustering_input)
  
  parse_dimension <- function(x, arg_name, allowed_units) {
    if (is.null(x)) {
      return(NULL)
    }
    
    if (length(x) != 2) {
      stop(sprintf(
        "'%s' must be NULL or a vector of length 2, e.g. c(14, \"in\").",
        arg_name
      ))
    }
    
    value <- suppressWarnings(as.numeric(x[1]))
    unit  <- as.character(x[2])
    
    if (is.na(value) || length(value) != 1 || value <= 0) {
      stop(sprintf(
        "The first element of '%s' must be a single positive number.",
        arg_name
      ))
    }
    
    if (!unit %in% allowed_units) {
      stop(sprintf(
        "The second element of '%s' must be one of: %s.",
        arg_name,
        paste(allowed_units, collapse = ", ")
      ))
    }
    
    list(value = value, unit = unit)
  }
  
  convert_to_inches <- function(value, unit) {
    if (unit == "in") return(value)
    if (unit == "cm") return(value / 2.54)
    if (unit == "mm") return(value / 25.4)
    if (unit == "px") return(value / 72)
    stop(sprintf("Unsupported unit for SVG conversion: %s", unit))
  }
  
  convert_png_dimension <- function(value, unit) {
    if (unit == "px") {
      return(list(value = round(value), units = "px"))
    }
    if (unit == "in") {
      return(list(value = value, units = "in"))
    }
    if (unit == "cm") {
      return(list(value = value / 2.54, units = "in"))
    }
    if (unit == "mm") {
      return(list(value = value / 25.4, units = "in"))
    }
    stop(sprintf("Unsupported unit for PNG conversion: %s", unit))
  }
  
  get_default_branch_colors <- function(k) {
    base_colors <- c(
      "#D73027", "#4575B4", "#1A9850", "#984EA3",
      "#FF8C00", "#4DAF4A", "#A65628", "#F781BF",
      "#999999", "#66C2A5", "#FC8D62", "#8DA0CB"
    )
    
    if (k <= length(base_colors)) {
      return(base_colors[seq_len(k)])
    }
    
    grDevices::colorRampPalette(base_colors)(k)
  }
  
  svg_width_parsed  <- parse_dimension(svg_width,  "svg_width",  valid_svg_units)
  svg_height_parsed <- parse_dimension(svg_height, "svg_height", valid_svg_units)
  png_width_parsed  <- parse_dimension(png_width,  "png_width",  valid_png_units)
  png_height_parsed <- parse_dimension(png_height, "png_height", valid_png_units)
  
  if (!distance_method %in% available_distance_methods) {
    stop(
      "Invalid 'distance_method'. Available options: ",
      paste(available_distance_methods, collapse = ", ")
    )
  }
  
  if (!clustering_method %in% available_clustering_methods) {
    stop(
      "Invalid 'clustering_method'. Available options: ",
      paste(available_clustering_methods, collapse = ", ")
    )
  }
  
  if (!all(label_fields_rows %in% available_label_fields)) {
    stop(
      "Invalid values in 'label_fields_rows'. Available options: ",
      paste(available_label_fields, collapse = ", ")
    )
  }
  
  if (!all(label_fields_cols %in% available_label_fields)) {
    stop(
      "Invalid values in 'label_fields_cols'. Available options: ",
      paste(available_label_fields, collapse = ", ")
    )
  }
  
  if (!is.logical(show_rg_values) || length(show_rg_values) != 1 || is.na(show_rg_values)) {
    stop("'show_rg_values' must be a single TRUE/FALSE value.")
  }
  
  if (!is.logical(show_significance_stars) || length(show_significance_stars) != 1 || is.na(show_significance_stars)) {
    stop("'show_significance_stars' must be a single TRUE/FALSE value.")
  }
  
  if (!is.logical(color_row_dend_branches) || length(color_row_dend_branches) != 1 || is.na(color_row_dend_branches)) {
    stop("'color_row_dend_branches' must be a single TRUE/FALSE value.")
  }
  
  if (!is.logical(color_col_dend_branches) || length(color_col_dend_branches) != 1 || is.na(color_col_dend_branches)) {
    stop("'color_col_dend_branches' must be a single TRUE/FALSE value.")
  }
  
  if (!is.logical(use_row_cluster_gaps) || length(use_row_cluster_gaps) != 1 || is.na(use_row_cluster_gaps)) {
    stop("'use_row_cluster_gaps' must be a single TRUE/FALSE value.")
  }
  
  if (!is.logical(use_col_cluster_gaps) || length(use_col_cluster_gaps) != 1 || is.na(use_col_cluster_gaps)) {
    stop("'use_col_cluster_gaps' must be a single TRUE/FALSE value.")
  }
  
  if (!is.character(row_gap_color) || length(row_gap_color) != 1 || is.na(row_gap_color)) {
    stop("'row_gap_color' must be a single character color value.")
  }
  
  if (!is.character(col_gap_color) || length(col_gap_color) != 1 || is.na(col_gap_color)) {
    stop("'col_gap_color' must be a single character color value.")
  }
  
  if (!is.numeric(row_gap_width_mm) || length(row_gap_width_mm) != 1 || is.na(row_gap_width_mm) || row_gap_width_mm < 0) {
    stop("'row_gap_width_mm' must be a single non-negative numeric value.")
  }
  
  if (!is.numeric(col_gap_width_mm) || length(col_gap_width_mm) != 1 || is.na(col_gap_width_mm) || col_gap_width_mm < 0) {
    stop("'col_gap_width_mm' must be a single non-negative numeric value.")
  }
  
  if (!is.null(row_dend_k)) {
    if (!is.numeric(row_dend_k) || length(row_dend_k) != 1 || is.na(row_dend_k) || row_dend_k < 1) {
      stop("'row_dend_k' must be NULL or a single positive integer.")
    }
    row_dend_k <- as.integer(row_dend_k)
  }
  
  if (!is.null(col_dend_k)) {
    if (!is.numeric(col_dend_k) || length(col_dend_k) != 1 || is.na(col_dend_k) || col_dend_k < 1) {
      stop("'col_dend_k' must be NULL or a single positive integer.")
    }
    col_dend_k <- as.integer(col_dend_k)
  }
  
  if (!is.null(row_dend_colors) && !is.character(row_dend_colors)) {
    stop("'row_dend_colors' must be NULL or a character vector of colors.")
  }
  
  if (!is.null(col_dend_colors) && !is.character(col_dend_colors)) {
    stop("'col_dend_colors' must be NULL or a character vector of colors.")
  }
  
  if (isTRUE(color_row_dend_branches) && !isTRUE(cluster)) {
    stop("'color_row_dend_branches = TRUE' requires 'cluster = TRUE'.")
  }
  
  if (isTRUE(color_col_dend_branches) && !isTRUE(cluster)) {
    stop("'color_col_dend_branches = TRUE' requires 'cluster = TRUE'.")
  }
  
  if (isTRUE(use_row_cluster_gaps) && !isTRUE(cluster)) {
    stop("'use_row_cluster_gaps = TRUE' requires 'cluster = TRUE'.")
  }
  
  if (isTRUE(use_col_cluster_gaps) && !isTRUE(cluster)) {
    stop("'use_col_cluster_gaps = TRUE' requires 'cluster = TRUE'.")
  }
  
  if (isTRUE(color_row_dend_branches) && is.null(row_dend_k)) {
    stop("When 'color_row_dend_branches = TRUE', 'row_dend_k' must be provided.")
  }
  
  if (isTRUE(color_col_dend_branches) && is.null(col_dend_k)) {
    stop("When 'color_col_dend_branches = TRUE', 'col_dend_k' must be provided.")
  }
  
  if (isTRUE(use_row_cluster_gaps) && is.null(row_dend_k)) {
    stop("When 'use_row_cluster_gaps = TRUE', 'row_dend_k' must be provided.")
  }
  
  if (isTRUE(use_col_cluster_gaps) && is.null(col_dend_k)) {
    stop("When 'use_col_cluster_gaps = TRUE', 'col_dend_k' must be provided.")
  }
  
  if (!is.null(output_svg) && !is.character(output_svg)) {
    stop("'output_svg' must be NULL or a character path.")
  }
  
  if (!is.null(output_png) && !is.character(output_png)) {
    stop("'output_png' must be NULL or a character path.")
  }
  
  if (!is.numeric(png_res) || length(png_res) != 1 || png_res <= 0) {
    stop("'png_res' must be a single positive numeric value.")
  }
  
  selected_traits_local <- unique(as.character(selected_traits))
  
  if (length(selected_traits_local) == 0) {
    stop("'selected_traits' must contain at least one trait ID.")
  }
  
  if (any(is.na(selected_traits_local)) || any(selected_traits_local == "")) {
    stop("'selected_traits' contains NA or empty values.")
  }
  
  force(selected_traits_local)
  
  log_step("Argument validation and package loading")
  
  # =========================================
  # 2. Pairwise data
  # =========================================
  pairwise_rg <- ldsc_results_annotated %>%
    filter(
      .data$p1_id %in% .env$selected_traits_local,
      .data$p2_id %in% .env$selected_traits_local
    ) %>%
    filter(!is.na(.data$summary_rg)) %>%
    select(.data$p1_id, .data$p2_id, .data$summary_rg, .data$summary_rg_p)
  
  pairwise_unique <- pairwise_rg %>%
    rowwise() %>%
    mutate(
      trait_min = sort(c(.data$p1_id, .data$p2_id))[1],
      trait_max = sort(c(.data$p1_id, .data$p2_id))[2]
    ) %>%
    ungroup() %>%
    group_by(.data$trait_min, .data$trait_max) %>%
    summarise(
      summary_rg = first(na.omit(.data$summary_rg)),
      summary_rg_p = first(na.omit(.data$summary_rg_p)),
      .groups = "drop"
    )
  
  pairwise_sym <- bind_rows(
    pairwise_unique %>%
      transmute(
        p1_id = .data$trait_min,
        p2_id = .data$trait_max,
        summary_rg = .data$summary_rg,
        summary_rg_p = .data$summary_rg_p
      ),
    pairwise_unique %>%
      transmute(
        p1_id = .data$trait_max,
        p2_id = .data$trait_min,
        summary_rg = .data$summary_rg,
        summary_rg_p = .data$summary_rg_p
      )
  )
  
  full_grid <- expand.grid(
    p1_id = selected_traits_local,
    p2_id = selected_traits_local,
    stringsAsFactors = FALSE
  ) %>%
    as_tibble()
  
  pairwise_rg_full <- full_grid %>%
    left_join(pairwise_sym, by = c("p1_id", "p2_id")) %>%
    mutate(
      summary_rg = ifelse(.data$p1_id == .data$p2_id, 1, .data$summary_rg),
      summary_rg_p = ifelse(.data$p1_id == .data$p2_id, NA, .data$summary_rg_p)
    )
  
  log_step("Preparing pairwise rg data and closing the full grid")
  
  # =========================================
  # 3. Build matrix + repair NA/NaN/Inf
  # =========================================
  rg_matrix <- pairwise_rg_full %>%
    pivot_wider(
      id_cols = .data$p1_id,
      names_from = .data$p2_id,
      values_from = .data$summary_rg
    ) %>%
    column_to_rownames("p1_id") %>%
    as.matrix()
  
  p_matrix <- pairwise_rg_full %>%
    pivot_wider(
      id_cols = .data$p1_id,
      names_from = .data$p2_id,
      values_from = .data$summary_rg_p
    ) %>%
    column_to_rownames("p1_id") %>%
    as.matrix()
  
  n_bad_rg <- sum(!is.finite(rg_matrix))
  n_bad_p  <- sum(!is.finite(p_matrix))
  
  rg_matrix[!is.finite(rg_matrix)] <- 0
  p_matrix[!is.finite(p_matrix)] <- 1
  
  log_step("Building rg and p-value matrices")
  if (isTRUE(verbose)) {
    cat(sprintf("           NA/NaN/Inf replaced: rg=%d, p=%d\n", n_bad_rg, n_bad_p))
  }
  
  # =========================================
  # 4. Matrix for display and clustering
  # =========================================
  heatmap_matrix <- if (absolute_correlation) abs(rg_matrix) else rg_matrix
  
  if (clustering_input == "signed_rg") {
    clustering_matrix <- rg_matrix
  } else if (clustering_input == "absolute_rg") {
    clustering_matrix <- abs(rg_matrix)
  } else if (clustering_input == "displayed_matrix") {
    clustering_matrix <- heatmap_matrix
  } else {
    stop("Internal error: unsupported 'clustering_input'.")
  }
  
  clustering_matrix[!is.finite(clustering_matrix)] <- 0
  
  row_split <- NULL
  col_split <- NULL
  
  if (cluster) {
    row_dist <- dist(clustering_matrix, method = distance_method)
    col_dist <- dist(t(clustering_matrix), method = distance_method)
    
    row_hc <- hclust(row_dist, method = clustering_method)
    col_hc <- hclust(col_dist, method = clustering_method)
    
    row_dend <- as.dendrogram(row_hc)
    col_dend <- as.dendrogram(col_hc)
    
    if (isTRUE(color_row_dend_branches)) {
      if (is.null(row_dend_colors)) {
        row_dend_colors <- get_default_branch_colors(row_dend_k)
      }
      
      row_dend <- dendextend::color_branches(
        dend = row_dend,
        k = row_dend_k,
        col = row_dend_colors
      )
    }
    
    if (isTRUE(color_col_dend_branches)) {
      if (is.null(col_dend_colors)) {
        col_dend_colors <- get_default_branch_colors(col_dend_k)
      }
      
      col_dend <- dendextend::color_branches(
        dend = col_dend,
        k = col_dend_k,
        col = col_dend_colors
      )
    }
    
    row_order_ids <- labels(row_dend)[order.dendrogram(row_dend)]
    col_order_ids <- labels(col_dend)[order.dendrogram(col_dend)]
    
    if (isTRUE(use_row_cluster_gaps)) {
      row_clusters <- cutree(row_hc, k = row_dend_k)
      row_split <- factor(
        row_clusters[rownames(heatmap_matrix)],
        levels = unique(row_clusters[row_order_ids])
      )
    }
    
    if (isTRUE(use_col_cluster_gaps)) {
      col_clusters <- cutree(col_hc, k = col_dend_k)
      col_split <- factor(
        col_clusters[colnames(heatmap_matrix)],
        levels = unique(col_clusters[col_order_ids])
      )
    }
  } else {
    row_hc <- FALSE
    col_hc <- FALSE
    row_dend <- FALSE
    col_dend <- FALSE
    row_order_ids <- rownames(rg_matrix)
    col_order_ids <- colnames(rg_matrix)
  }
  
  log_step("Clustering and setting trait order")
  
  if (isTRUE(verbose)) {
    cat("           clustering_input:", clustering_input, "\n")
    cat("           distance_method:", distance_method, "\n")
    cat("           clustering_method:", clustering_method, "\n")
    
    if (cluster) {
      cat("           row order:\n")
      print(row_order_ids)
      cat("           column order:\n")
      print(col_order_ids)
      
      if (isTRUE(color_row_dend_branches)) {
        cat("           row dendrogram branches colored with k =", row_dend_k, "\n")
      }
      
      if (isTRUE(color_col_dend_branches)) {
        cat("           column dendrogram branches colored with k =", col_dend_k, "\n")
      }
      
      if (isTRUE(use_row_cluster_gaps)) {
        cat("           row cluster gaps enabled with k =", row_dend_k, "\n")
      }
      
      if (isTRUE(use_col_cluster_gaps)) {
        cat("           column cluster gaps enabled with k =", col_dend_k, "\n")
      }
    }
  }
  
  # =========================================
  # 5. Metadata for labels
  # =========================================
  metadata_base <- metadata_ieu_EURsampleSize10000 %>%
    select(any_of(c(
      "id",
      "trait",
      "category",
      "subcategory",
      "sample_size",
      "coverage",
      "year",
      "author",
      "sex",
      "population",
      "unit",
      "ncase",
      "ncontrol",
      "consortium",
      "study_design",
      "pmid",
      "doi",
      "ontology",
      "group_name",
      "build",
      "covariates",
      "priority",
      "mr",
      "qc_prior_to_upload",
      "note",
      "sd"
    ))) %>%
    distinct()
  
  ldsc_label_base <- ldsc_results_annotated %>%
    filter(!is.na(.data$p1_id)) %>%
    select(
      id = .data$p1_id,
      nsnp = .data$p1_read_snps
    ) %>%
    distinct() %>%
    group_by(.data$id) %>%
    summarise(
      nsnp = first(na.omit(.data$nsnp)),
      .groups = "drop"
    )
  
  label_metadata <- metadata_base %>%
    left_join(ldsc_label_base, by = "id") %>%
    select(any_of(available_label_fields)) %>%
    distinct()
  
  format_label_value <- function(field_name, value) {
    if (is.na(value)) {
      return("NA")
    }
    
    if (field_name == "id") return(as.character(value))
    if (field_name == "trait") return(as.character(value))
    if (field_name == "category") return(paste0("category: ", value))
    if (field_name == "subcategory") return(as.character(value))
    if (field_name == "nsnp") return(paste0("nsnp: ", value))
    if (field_name == "sample_size") return(paste0("n: ", value))
    if (field_name == "coverage") return(paste0("coverage: ", value))
    if (field_name == "year") return(paste0("year: ", value))
    if (field_name == "author") return(paste0("author: ", value))
    if (field_name == "sex") return(paste0("sex: ", value))
    if (field_name == "population") return(paste0("population: ", value))
    if (field_name == "unit") return(paste0("unit: ", value))
    if (field_name == "ncase") return(paste0("ncase: ", value))
    if (field_name == "ncontrol") return(paste0("ncontrol: ", value))
    if (field_name == "consortium") return(paste0("consortium: ", value))
    if (field_name == "study_design") return(paste0("study design: ", value))
    if (field_name == "pmid") return(paste0("PMID: ", value))
    if (field_name == "doi") return(paste0("DOI: ", value))
    if (field_name == "ontology") return(paste0("ontology: ", value))
    if (field_name == "group_name") return(paste0("group name: ", value))
    if (field_name == "build") return(paste0("build: ", value))
    if (field_name == "covariates") return(paste0("covariates: ", value))
    if (field_name == "priority") return(paste0("priority: ", value))
    if (field_name == "mr") return(paste0("mr: ", value))
    if (field_name == "qc_prior_to_upload") return(paste0("qc prior to upload: ", value))
    if (field_name == "note") return(paste0("note: ", value))
    if (field_name == "sd") return(paste0("sd: ", value))
    
    paste0(field_name, ": ", value)
  }
  
  build_label_vector <- function(label_fields, label_metadata, matrix_ids) {
    label_vector <- setNames(
      lapply(seq_len(nrow(label_metadata)), function(i) {
        row_i <- label_metadata[i, , drop = FALSE]
        
        pieces <- sapply(label_fields, function(field) {
          format_label_value(field, row_i[[field]][1])
        }, USE.NAMES = FALSE)
        
        paste(pieces, collapse = label_separator)
      }),
      label_metadata$id
    )
    
    label_vector <- unlist(label_vector)
    missing_ids <- setdiff(matrix_ids, names(label_vector))
    label_vector[missing_ids] <- missing_ids
    label_vector
  }
  
  row_label_vector <- build_label_vector(
    label_fields = label_fields_rows,
    label_metadata = label_metadata,
    matrix_ids = rownames(rg_matrix)
  )
  
  col_label_vector <- build_label_vector(
    label_fields = label_fields_cols,
    label_metadata = label_metadata,
    matrix_ids = colnames(rg_matrix)
  )
  
  row_labels <- row_label_vector[rownames(rg_matrix)]
  col_labels <- col_label_vector[colnames(rg_matrix)]
  
  log_step("Preparing row and column labels")
  
  # =========================================
  # 6. Significance stars
  # =========================================
  get_stars <- function(p) {
    ifelse(
      is.na(p), "",
      ifelse(
        p < 0.001, "***",
        ifelse(
          p < 0.01, "**",
          ifelse(p < 0.05, "*", "")
        )
      )
    )
  }
  
  stars_matrix <- apply(p_matrix, c(1, 2), get_stars)
  
  log_step("Calculating significance stars")
  
  # =========================================
  # 7. Color scale + text layout
  # =========================================
  if (absolute_correlation) {
    col_fun <- colorRamp2(c(0, 1), c("white", absolute_high_color))
    legend_title <- "Absolute genetic\ncorrelation |rg|"
    legend_at <- c(0, 0.25, 0.5, 0.75, 1)
  } else {
    col_fun <- colorRamp2(c(-1, 0, 1), c(low_color, mid_color, high_color))
    legend_title <- "Genetic\ncorrelation (rg)"
    legend_at <- c(-1, -0.5, 0, 0.5, 1)
  }
  
  max_row_w <- max_text_width(
    row_labels,
    gp = gpar(fontsize = label_size)
  )
  
  max_col_h <- max_text_width(
    col_labels,
    gp = gpar(fontsize = label_size)
  )
  
  row_names_max_width <- min(
    max_row_w + unit(3, "mm"),
    unit(row_names_max_width_mm_cap, "mm")
  )
  
  column_names_max_height <- min(
    max_col_h + unit(4, "mm"),
    unit(column_names_max_height_mm_cap, "mm")
  )
  
  get_text_color <- function(v) {
    text_color <- "black"
    
    if (!is.null(text_white_high_threshold) && !is.na(v) && v >= text_white_high_threshold) {
      text_color <- "white"
    }
    
    if (!is.null(text_white_low_threshold) && !is.na(v) && v <= text_white_low_threshold) {
      text_color <- "white"
    }
    
    text_color
  }
  
  cell_fun <- function(j, i, x, y, width, height, fill) {
    val <- rg_matrix[i, j]
    txt_col <- get_text_color(val)
    
    if (isTRUE(show_rg_values) && !is.na(val)) {
      grid.text(
        sprintf("%.2f", val),
        x = x,
        y = y + unit(rg_text_y_offset_mm, "mm"),
        gp = gpar(
          fontsize = rg_size,
          fontface = "bold",
          col = txt_col
        )
      )
    }
    
    if (isTRUE(show_significance_stars) && stars_matrix[i, j] != "") {
      grid.text(
        stars_matrix[i, j],
        x = x,
        y = y + unit(star_text_y_offset_mm, "mm"),
        gp = gpar(
          fontsize = star_size,
          fontface = "bold",
          col = txt_col
        )
      )
    }
  }
  
  log_step("Preparing colors, legend, and text layout")
  
  # =========================================
  # 8. Heatmap object + draw / save
  # =========================================
  ht <- Heatmap(
    heatmap_matrix,
    name = "rg",
    col = col_fun,
    
    rect_gp = gpar(col = "white", lwd = 0.8),
    
    cluster_rows = row_dend,
    cluster_columns = col_dend,
    show_row_dend = isTRUE(cluster),
    show_column_dend = isTRUE(cluster),
    
    row_split = row_split,
    column_split = col_split,
    row_gap = unit(row_gap_width_mm, "mm"),
    column_gap = unit(col_gap_width_mm, "mm"),
    
    row_dend_width = unit(row_dend_width_mm, "mm"),
    column_dend_height = unit(column_dend_height_mm, "mm"),
    
    row_names_max_width = row_names_max_width,
    column_names_max_height = column_names_max_height,
    
    row_labels = row_labels,
    column_labels = col_labels,
    
    row_names_side = row_labels_side,
    column_names_side = col_labels_side,
    
    row_names_gp = gpar(fontsize = label_size),
    column_names_gp = gpar(fontsize = label_size),
    column_names_rot = 90,
    
    cell_fun = cell_fun,
    
    heatmap_legend_param = list(
      title = legend_title,
      at = legend_at,
      legend_height = unit(55, "mm"),
      legend_width = unit(8, "mm"),
      title_gp = gpar(fontsize = 9, fontface = "bold"),
      labels_gp = gpar(fontsize = 8)
    ),
    
    column_title = title,
    column_title_gp = gpar(fontsize = 12, fontface = "bold")
  )
  
  draw_heatmap <- function() {
    draw(
      ht,
      padding = unit(
        c(
          top_padding_mm,
          right_padding_mm + legend_gap_mm,
          bottom_padding_mm,
          left_padding_mm
        ),
        "mm"
      ),
      heatmap_legend_side = "right"
    )
  }
  
  if (!is.null(output_svg)) {
    svg_dir <- dirname(output_svg)
    if (!dir.exists(svg_dir)) {
      dir.create(svg_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    svg_args <- list(filename = output_svg)
    
    if (!is.null(svg_width_parsed)) {
      svg_args$width <- convert_to_inches(svg_width_parsed$value, svg_width_parsed$unit)
    }
    
    if (!is.null(svg_height_parsed)) {
      svg_args$height <- convert_to_inches(svg_height_parsed$value, svg_height_parsed$unit)
    }
    
    do.call(svg, svg_args)
    draw_heatmap()
    dev.off()
    
    if (file.exists(output_svg)) {
      svg_size_mb <- file.info(output_svg)$size / (1024^2)
      cat(sprintf("Saved SVG: %s (%.2f MB)\n", output_svg, svg_size_mb))
    }
  }
  
  if (!is.null(output_png)) {
    png_dir <- dirname(output_png)
    if (!dir.exists(png_dir)) {
      dir.create(png_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    png_args <- list(
      filename = output_png,
      res = png_res
    )
    
    if (is.null(png_width_parsed) && is.null(png_height_parsed)) {
      png_args$width <- 10
      png_args$height <- 10
      png_args$units <- "in"
    } else if (!is.null(png_width_parsed) && is.null(png_height_parsed)) {
      width_conv <- convert_png_dimension(png_width_parsed$value, png_width_parsed$unit)
      png_args$width <- width_conv$value
      png_args$height <- 10
      png_args$units <- width_conv$units
    } else if (is.null(png_width_parsed) && !is.null(png_height_parsed)) {
      height_conv <- convert_png_dimension(png_height_parsed$value, png_height_parsed$unit)
      png_args$width <- 10
      png_args$height <- height_conv$value
      png_args$units <- height_conv$units
    } else {
      width_conv  <- convert_png_dimension(png_width_parsed$value, png_width_parsed$unit)
      height_conv <- convert_png_dimension(png_height_parsed$value, png_height_parsed$unit)
      
      if (width_conv$units != height_conv$units) {
        stop("For PNG, width and height must resolve to the same unit type.")
      }
      
      png_args$width <- width_conv$value
      png_args$height <- height_conv$value
      png_args$units <- width_conv$units
    }
    
    do.call(png, png_args)
    draw_heatmap()
    dev.off()
    
    if (file.exists(output_png)) {
      png_size_mb <- file.info(output_png)$size / (1024^2)
      cat(sprintf("Saved PNG: %s (%.2f MB)\n", output_png, png_size_mb))
    }
  }
  
  if (is.null(output_svg) && is.null(output_png)) {
    draw_heatmap()
  }
  
  log_step("Rendering heatmap and optional saving to SVG/PNG")
  
  # =========================================
  # Return
  # =========================================
  if (isTRUE(return_data)) {
    return(list(
      heatmap = ht,
      rg_matrix = rg_matrix,
      p_matrix = p_matrix,
      heatmap_matrix = heatmap_matrix,
      clustering_matrix = clustering_matrix,
      row_labels = row_labels,
      col_labels = col_labels,
      label_metadata = label_metadata,
      stars_matrix = stars_matrix,
      row_hclust = row_hc,
      col_hclust = col_hc,
      row_dendrogram = row_dend,
      col_dendrogram = col_dend,
      row_split = row_split,
      col_split = col_split,
      row_order = row_order_ids,
      col_order = col_order_ids,
      row_gap_color = row_gap_color,
      col_gap_color = col_gap_color,
      selected_traits = selected_traits_local
    ))
  } else {
    invisible(ht)
  }
}