create_rg_heatmap <- function(
    selected_traits,
    ldsc_results_annotated,
    metadata_ieu_EURsampleSize10000,
    
    title = "Genetic correlation heatmap (LDSC)",
    
    label_size = 8.5,
    rg_size = 7.5,
    star_size = 9,
    
    # 🔥 NOWE: offset tekstu
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
    text_white_low_threshold = NULL
) {
  
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(grid)
  
  # =========================================
  # 1. Walidacja
  # =========================================
  available_distance_methods <- c(
    "euclidean", "maximum", "manhattan",
    "canberra", "binary", "minkowski"
  )
  
  available_clustering_methods <- c(
    "complete", "single", "average", "mcquitty",
    "median", "centroid", "ward.D", "ward.D2"
  )
  
  if (!distance_method %in% available_distance_methods) {
    stop("Wrong distance_method")
  }
  
  if (!clustering_method %in% available_clustering_methods) {
    stop("Wrong clustering_method")
  }
  
  # =========================================
  # 2. Dane
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
      transmute(p1_id = trait_min, p2_id = trait_max, summary_rg, summary_rg_p),
    pairwise_unique %>%
      transmute(p1_id = trait_max, p2_id = trait_min, summary_rg, summary_rg_p)
  )
  
  full_grid <- expand.grid(
    p1_id = selected_traits,
    p2_id = selected_traits,
    stringsAsFactors = FALSE
  ) %>% as_tibble()
  
  pairwise_rg_full <- full_grid %>%
    left_join(pairwise_sym, by = c("p1_id", "p2_id")) %>%
    mutate(
      summary_rg = ifelse(p1_id == p2_id, 1, summary_rg),
      summary_rg_p = ifelse(p1_id == p2_id, NA, summary_rg_p)
    )
  
  rg_matrix <- pairwise_rg_full %>%
    pivot_wider(id_cols = p1_id, names_from = p2_id, values_from = summary_rg) %>%
    column_to_rownames("p1_id") %>%
    as.matrix()
  
  p_matrix <- pairwise_rg_full %>%
    pivot_wider(id_cols = p1_id, names_from = p2_id, values_from = summary_rg_p) %>%
    column_to_rownames("p1_id") %>%
    as.matrix()
  
  # =========================================
  # 3. abs mode
  # =========================================
  heatmap_matrix <- if (absolute_correlation) abs(rg_matrix) else rg_matrix
  
  rg_for_clust <- heatmap_matrix
  rg_for_clust[is.na(rg_for_clust)] <- 0
  
  # =========================================
  # 4. clustering
  # =========================================
  if (cluster) {
    hc <- hclust(
      dist(rg_for_clust, method = distance_method),
      method = clustering_method
    )
    ord <- hc$labels[hc$order]
  } else {
    ord <- rownames(rg_matrix)
  }
  
  rg_matrix <- rg_matrix[ord, ord]
  heatmap_matrix <- heatmap_matrix[ord, ord]
  p_matrix <- p_matrix[ord, ord]
  
  # =========================================
  # 5. labels
  # =========================================
  label_vector <- metadata_ieu_EURsampleSize10000 %>%
    select(id, trait) %>%
    distinct() %>%
    mutate(label = paste0(id, " | ", trait)) %>%
    { setNames(.$label, .$id) }
  
  row_labels <- label_vector[rownames(rg_matrix)]
  col_labels <- label_vector[colnames(rg_matrix)]
  
  # =========================================
  # 6. stars
  # =========================================
  get_stars <- function(p) {
    ifelse(is.na(p), "",
           ifelse(p < 0.001, "***",
                  ifelse(p < 0.01, "**",
                         ifelse(p < 0.05, "*", ""))))
  }
  
  stars_matrix <- apply(p_matrix, c(1,2), get_stars)
  
  # =========================================
  # 7. colors
  # =========================================
  if (absolute_correlation) {
    col_fun <- colorRamp2(c(0,1), c("white", absolute_high_color))
  } else {
    col_fun <- colorRamp2(c(-1,0,1), c(low_color, mid_color, high_color))
  }
  
  # =========================================
  # 8. text color logic
  # =========================================
  get_text_color <- function(v) {
    col <- "black"
    if (!is.null(text_white_high_threshold) && v >= text_white_high_threshold) col <- "white"
    if (!is.null(text_white_low_threshold) && v <= text_white_low_threshold) col <- "white"
    col
  }
  
  # =========================================
  # 9. cell_fun (🔥 tu jest nowa kontrola)
  # =========================================
  cell_fun <- function(j, i, x, y, width, height, fill) {
    
    val <- rg_matrix[i, j]
    txt_col <- get_text_color(val)
    
    if (!is.na(val)) {
      grid.text(
        sprintf("%.2f", val),
        x = x,
        y = y + unit(rg_text_y_offset_mm, "mm"),
        gp = gpar(fontsize = rg_size, fontface = "bold", col = txt_col)
      )
    }
    
    if (stars_matrix[i, j] != "") {
      grid.text(
        stars_matrix[i, j],
        x = x,
        y = y + unit(star_text_y_offset_mm, "mm"),
        gp = gpar(fontsize = star_size, col = txt_col, fontface = "bold")
      )
    }
  }
  
  # =========================================
  # 10. heatmap
  # =========================================
  ht <- Heatmap(
    heatmap_matrix,
    col = col_fun,
    rect_gp = gpar(col = "white"),
    
    cluster_rows = cluster,
    cluster_columns = cluster,
    
    row_labels = row_labels,
    column_labels = col_labels,
    
    row_names_side = row_labels_side,
    column_names_side = col_labels_side,
    
    row_names_gp = gpar(fontsize = label_size),
    column_names_gp = gpar(fontsize = label_size),
    column_names_rot = 90,
    
    row_dend_width = unit(28, "mm"),
    column_dend_height = unit(24, "mm"),
    
    cell_fun = cell_fun,
    
    heatmap_legend_param = list(
      title = ifelse(absolute_correlation, "|rg|", "rg"),
      legend_height = unit(50, "mm")
    ),
    
    column_title = title
  )
  
  # =========================================
  # 11. draw
  # =========================================
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
      )
    )
  }
  
  # =========================================
  # 12. svg save
  # =========================================
  if (!is.null(output_svg)) {
    
    if (is.null(svg_width) && is.null(svg_height)) {
      svg(output_svg)
    } else {
      svg(output_svg, width = svg_width, height = svg_height)
    }
    
    draw_heatmap()
    dev.off()
    
  } else {
    draw_heatmap()
  }
  
  invisible(ht)
}


create_rg_heatmap(
  selected_traits = selected_traits,
  ldsc_results_annotated = ldsc_results_annotated,
  metadata_ieu_EURsampleSize10000 = metadata_ieu_EURsampleSize10000,
  absolute_correlation = TRUE,
  text_white_high_threshold = 0.9,
  text_white_low_threshold = -0.9,
  rg_text_y_offset_mm = 1.0,
  star_text_y_offset_mm = -2.5,
)
