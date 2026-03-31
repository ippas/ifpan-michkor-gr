create_rg_heatmap <- function(
    selected_traits,
    ldsc_results_annotated,
    metadata_ieu_EURsampleSize10000,
    
    title = "Genetic correlation heatmap (LDSC)",
    
    label_size = 8.5,
    rg_size = 7.5,
    star_size = 9,
    
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
    
    legend_gap_mm = 10,
    
    top_padding_mm = 8,
    right_padding_mm = 8,
    bottom_padding_mm = 8,
    left_padding_mm = 8,
    
    output_svg = NULL,
    svg_width = NULL,
    svg_height = NULL,
    
    text_white_high_threshold = NULL,
    text_white_low_threshold = NULL,
    
    row_dend_width_mm = 28,
    column_dend_height_mm = 24,
    
    row_names_max_width_mm_cap = 85,
    column_names_max_height_mm_cap = 70,
    
    label_fields_rows = c("id", "trait"),
    label_fields_cols = c("id", "trait"),
    label_separator = " | ",
    
    verbose = TRUE,
    return_data = FALSE
) {
  
  # =========================================
  # 0. Helper do komunikatów
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
  # 1. Pakiety + walidacja
  # =========================================
  suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(circlize)
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(grid)
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
  
  if (!distance_method %in% available_distance_methods) {
    stop(
      "Nieprawidłowy 'distance_method'. Dostępne opcje: ",
      paste(available_distance_methods, collapse = ", ")
    )
  }
  
  if (!clustering_method %in% available_clustering_methods) {
    stop(
      "Nieprawidłowy 'clustering_method'. Dostępne opcje: ",
      paste(available_clustering_methods, collapse = ", ")
    )
  }
  
  if (!all(label_fields_rows %in% available_label_fields)) {
    stop(
      "Nieprawidłowe wartości w 'label_fields_rows'. Dostępne opcje: ",
      paste(available_label_fields, collapse = ", ")
    )
  }
  
  if (!all(label_fields_cols %in% available_label_fields)) {
    stop(
      "Nieprawidłowe wartości w 'label_fields_cols'. Dostępne opcje: ",
      paste(available_label_fields, collapse = ", ")
    )
  }
  
  if (!is.null(output_svg) && !is.character(output_svg)) {
    stop("'output_svg' musi być NULL albo ścieżką tekstową.")
  }
  
  if (!is.null(svg_width) && (!is.numeric(svg_width) || length(svg_width) != 1)) {
    stop("'svg_width' musi być NULL albo pojedynczą wartością numeryczną.")
  }
  
  if (!is.null(svg_height) && (!is.numeric(svg_height) || length(svg_height) != 1)) {
    stop("'svg_height' musi być NULL albo pojedynczą wartością numeryczną.")
  }
  
  log_step("Walidacja argumentów i ładowanie pakietów")
  
  # =========================================
  # 2. Pairwise data
  # =========================================
  pairwise_rg <- ldsc_results_annotated %>%
    filter(p1_id %in% selected_traits, p2_id %in% selected_traits) %>%
    filter(!is.na(summary_rg)) %>%
    select(p1_id, p2_id, summary_rg, summary_rg_p)
  
  pairwise_unique <- pairwise_rg %>%
    rowwise() %>%
    mutate(
      trait_min = sort(c(p1_id, p2_id))[1],
      trait_max = sort(c(p1_id, p2_id))[2]
    ) %>%
    ungroup() %>%
    group_by(trait_min, trait_max) %>%
    summarise(
      summary_rg = first(na.omit(summary_rg)),
      summary_rg_p = first(na.omit(summary_rg_p)),
      .groups = "drop"
    )
  
  pairwise_sym <- bind_rows(
    pairwise_unique %>%
      transmute(
        p1_id = trait_min,
        p2_id = trait_max,
        summary_rg,
        summary_rg_p
      ),
    pairwise_unique %>%
      transmute(
        p1_id = trait_max,
        p2_id = trait_min,
        summary_rg,
        summary_rg_p
      )
  )
  
  full_grid <- expand.grid(
    p1_id = selected_traits,
    p2_id = selected_traits,
    stringsAsFactors = FALSE
  ) %>%
    as_tibble()
  
  pairwise_rg_full <- full_grid %>%
    left_join(pairwise_sym, by = c("p1_id", "p2_id")) %>%
    mutate(
      summary_rg = ifelse(p1_id == p2_id, 1, summary_rg),
      summary_rg_p = ifelse(p1_id == p2_id, NA, summary_rg_p)
    )
  
  log_step("Przygotowanie pairwise rg i domknięcie pełnej siatki")
  
  # =========================================
  # 3. Budowa macierzy
  # =========================================
  rg_matrix <- pairwise_rg_full %>%
    pivot_wider(
      id_cols = p1_id,
      names_from = p2_id,
      values_from = summary_rg
    ) %>%
    column_to_rownames("p1_id") %>%
    as.matrix()
  
  p_matrix <- pairwise_rg_full %>%
    pivot_wider(
      id_cols = p1_id,
      names_from = p2_id,
      values_from = summary_rg_p
    ) %>%
    column_to_rownames("p1_id") %>%
    as.matrix()
  
  log_step("Budowa macierzy rg i p-value")
  
  # =========================================
  # 4. Macierz do kolorów i klastrowania
  # =========================================
  heatmap_matrix <- if (absolute_correlation) abs(rg_matrix) else rg_matrix
  
  rg_for_clust <- heatmap_matrix
  rg_for_clust[is.na(rg_for_clust)] <- 0
  
  if (cluster) {
    hc <- hclust(
      dist(rg_for_clust, method = distance_method),
      method = clustering_method
    )
    ord <- hc$labels[hc$order]
  } else {
    ord <- rownames(rg_matrix)
  }
  
  rg_matrix <- rg_matrix[ord, ord, drop = FALSE]
  heatmap_matrix <- heatmap_matrix[ord, ord, drop = FALSE]
  p_matrix <- p_matrix[ord, ord, drop = FALSE]
  
  log_step("Klastrowanie i ustawienie kolejności traitów")
  
  # =========================================
  # 5. Metadane do etykiet
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
    filter(!is.na(p1_id)) %>%
    select(
      id = p1_id,
      nsnp = p1_read_snps
    ) %>%
    distinct() %>%
    group_by(id) %>%
    summarise(
      nsnp = first(na.omit(nsnp)),
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
    
    if (field_name == "id") {
      return(as.character(value))
    }
    
    if (field_name == "trait") {
      return(as.character(value))
    }
    
    if (field_name == "category") {
      return(paste0("category: ", value))
    }
    
    if (field_name == "subcategory") {
      return(as.character(value))
    }
    
    if (field_name == "nsnp") {
      return(paste0("nsnp: ", value))
    }
    
    if (field_name == "sample_size") {
      return(paste0("n: ", value))
    }
    
    if (field_name == "coverage") {
      return(paste0("coverage: ", value))
    }
    
    if (field_name == "year") {
      return(paste0("year: ", value))
    }
    
    if (field_name == "author") {
      return(paste0("author: ", value))
    }
    
    if (field_name == "sex") {
      return(paste0("sex: ", value))
    }
    
    if (field_name == "population") {
      return(paste0("population: ", value))
    }
    
    if (field_name == "unit") {
      return(paste0("unit: ", value))
    }
    
    if (field_name == "ncase") {
      return(paste0("ncase: ", value))
    }
    
    if (field_name == "ncontrol") {
      return(paste0("ncontrol: ", value))
    }
    
    if (field_name == "consortium") {
      return(paste0("consortium: ", value))
    }
    
    if (field_name == "study_design") {
      return(paste0("study design: ", value))
    }
    
    if (field_name == "pmid") {
      return(paste0("PMID: ", value))
    }
    
    if (field_name == "doi") {
      return(paste0("DOI: ", value))
    }
    
    if (field_name == "ontology") {
      return(paste0("ontology: ", value))
    }
    
    if (field_name == "group_name") {
      return(paste0("group name: ", value))
    }
    
    if (field_name == "build") {
      return(paste0("build: ", value))
    }
    
    if (field_name == "covariates") {
      return(paste0("covariates: ", value))
    }
    
    if (field_name == "priority") {
      return(paste0("priority: ", value))
    }
    
    if (field_name == "mr") {
      return(paste0("mr: ", value))
    }
    
    if (field_name == "qc_prior_to_upload") {
      return(paste0("qc prior to upload: ", value))
    }
    
    if (field_name == "note") {
      return(paste0("note: ", value))
    }
    
    if (field_name == "sd") {
      return(paste0("sd: ", value))
    }
    
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
  
  log_step("Przygotowanie etykiet wierszy i kolumn")
  
  # =========================================
  # 6. Gwiazdki istotności
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
  
  log_step("Obliczenie gwiazdek istotności")
  
  # =========================================
  # 7. Skala kolorów + layout tekstu
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
    
    if (!is.na(val)) {
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
    
    if (stars_matrix[i, j] != "") {
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
  
  log_step("Przygotowanie kolorów, legendy i layoutu tekstu")
  
  # =========================================
  # 8. Heatmap object + rysowanie / zapis
  # =========================================
  ht <- Heatmap(
    heatmap_matrix,
    name = "rg",
    col = col_fun,
    
    rect_gp = gpar(col = "white", lwd = 0.8),
    
    cluster_rows = cluster,
    cluster_columns = cluster,
    show_row_dend = cluster,
    show_column_dend = cluster,
    
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
    
    if (is.null(svg_width) && is.null(svg_height)) {
      svg(filename = output_svg)
    } else if (!is.null(svg_width) && is.null(svg_height)) {
      svg(filename = output_svg, width = svg_width)
    } else if (is.null(svg_width) && !is.null(svg_height)) {
      svg(filename = output_svg, height = svg_height)
    } else {
      svg(filename = output_svg, width = svg_width, height = svg_height)
    }
    
    draw_heatmap()
    dev.off()
  } else {
    draw_heatmap()
  }
  
  log_step("Renderowanie heatmapy i opcjonalny zapis do SVG")
  
  # =========================================
  # Return
  # =========================================
  if (isTRUE(return_data)) {
    return(list(
      heatmap = ht,
      rg_matrix = rg_matrix,
      p_matrix = p_matrix,
      heatmap_matrix = heatmap_matrix,
      row_labels = row_labels,
      col_labels = col_labels,
      label_metadata = label_metadata,
      stars_matrix = stars_matrix,
      order = ord
    ))
  } else {
    invisible(ht)
  }
}

# create_rg_heatmap(
#   selected_traits = selected_traits,
#   ldsc_results_annotated = ldsc_results_annotated,
#   metadata_ieu_EURsampleSize10000 = metadata_ieu_EURsampleSize10000,
#   absolute_correlation = TRUE,
#   text_white_high_threshold = 0.9,
#   text_white_low_threshold = -0.9,
#   rg_text_y_offset_mm = 1.0,
#   star_text_y_offset_mm = -2.5,
#   label_fields_rows = c("id", "trait", "subcategory", "nsnp", "sample_size"),
#   label_fields_cols = c("id", "trait"),
#   
#   row_dend_width_mm = 28,
#   column_dend_height_mm = 24,
#   
#   row_names_max_width_mm_cap = 185,
#   column_names_max_height_mm_cap = 170,
# )
