create_rg_heatmap <- function(
    selected_traits,
    ldsc_results_annotated,
    metadata_ieu_EURsampleSize10000,
    
    title = "Genetic correlation heatmap (LDSC)",
    
    label_size = 8.5,
    rg_size = 7.5,
    star_size = 9,
    
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
    left_padding_mm = 8
) {
  
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(grid)
  
  # =========================================
  # 1. Walidacja metod klastrowania
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
  
  # =========================================
  # 2. Przygotowanie danych pairwise rg
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
  
  # =========================================
  # 3. Budowa macierzy rg i p
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
  
  # =========================================
  # 4. Macierz do kolorów i klastrowania
  # =========================================
  if (absolute_correlation) {
    heatmap_matrix <- abs(rg_matrix)
  } else {
    heatmap_matrix <- rg_matrix
  }
  
  rg_for_clust <- heatmap_matrix
  rg_for_clust[is.na(rg_for_clust)] <- 0
  
  # =========================================
  # 5. Klastrowanie i kolejność
  # =========================================
  if (cluster) {
    hc <- hclust(
      dist(rg_for_clust, method = distance_method),
      method = clustering_method
    )
    trait_order <- hc$labels[hc$order]
  } else {
    trait_order <- rownames(rg_matrix)
  }
  
  rg_matrix <- rg_matrix[trait_order, trait_order, drop = FALSE]
  heatmap_matrix <- heatmap_matrix[trait_order, trait_order, drop = FALSE]
  p_matrix  <- p_matrix[trait_order, trait_order, drop = FALSE]
  
  # =========================================
  # 6. Etykiety
  # =========================================
  label_vector <- metadata_ieu_EURsampleSize10000 %>%
    select(id, trait) %>%
    distinct() %>%
    mutate(label = paste0(id, " | ", trait)) %>%
    { setNames(.$label, .$id) }
  
  missing_ids <- setdiff(rownames(rg_matrix), names(label_vector))
  label_vector[missing_ids] <- missing_ids
  
  row_labels <- label_vector[rownames(rg_matrix)]
  col_labels <- label_vector[colnames(rg_matrix)]
  
  # =========================================
  # 7. Gwiazdki istotności
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
  
  # =========================================
  # 8. Kolory
  # =========================================
  if (absolute_correlation) {
    col_fun <- colorRamp2(
      c(0, 1),
      c("#FFFFFF", absolute_high_color)
    )
    
    legend_title <- "Absolute genetic\ncorrelation |rg|"
    legend_at <- c(0, 0.25, 0.5, 0.75, 1)
  } else {
    col_fun <- colorRamp2(
      c(-1, 0, 1),
      c(low_color, mid_color, high_color)
    )
    
    legend_title <- "Genetic\ncorrelation (rg)"
    legend_at <- c(-1, -0.5, 0, 0.5, 1)
  }
  
  # =========================================
  # 9. Dynamiczne rozmiary layoutu
  # =========================================
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
    unit(85, "mm")
  )
  
  column_names_max_height <- min(
    max_col_h + unit(4, "mm"),
    unit(70, "mm")
  )
  
  # =========================================
  # 10. Tekst w komórkach
  # =========================================
  cell_fun <- function(j, i, x, y, width, height, fill) {
    
    if (!is.na(rg_matrix[i, j])) {
      grid.text(
        sprintf("%.2f", rg_matrix[i, j]),
        x = x,
        y = y + unit(1.0, "mm"),
        gp = gpar(
          fontsize = rg_size,
          fontface = "bold"
        )
      )
    }
    
    if (stars_matrix[i, j] != "") {
      grid.text(
        stars_matrix[i, j],
        x = x,
        y = y - unit(1.8, "mm"),
        gp = gpar(
          fontsize = star_size,
          col = "#1a1a1a",
          fontface = "bold"
        )
      )
    }
  }
  
  # =========================================
  # 11. Heatmap object
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
    
    row_dend_width = unit(28, "mm"),
    column_dend_height = unit(24, "mm"),
    
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
  
  # =========================================
  # 12. Rysowanie
  # =========================================
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
  
  invisible(ht)
}

create_rg_heatmap(
  selected_traits = selected_traits,
  ldsc_results_annotated = ldsc_results_annotated,
  metadata_ieu_EURsampleSize10000 = metadata_ieu_EURsampleSize10000,
  distance_method = "manhattan",
  clustering_method = "average"
)

create_rg_heatmap(
  selected_traits = selected_traits,
  ldsc_results_annotated = ldsc_results_annotated,
  metadata_ieu_EURsampleSize10000 = metadata_ieu_EURsampleSize10000,
  absolute_correlation = TRUE
)
