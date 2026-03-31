create_rg_heatmap <- function(
    selected_traits,
    ldsc_results_annotated,
    metadata_ieu_EURsampleSize10000,
    
    title = "Genetic correlation heatmap (LDSC)",
    
    label_size = 8.8,
    rg_size = 7.8,
    star_size = 9.5,
    
    low_color  = "#2166AC",
    mid_color  = "#F7F7F7",
    high_color = "#B2182B",
    
    row_labels_side = "left",
    col_labels_side = "top",
    
    cluster = TRUE,
    
    legend_gap_mm = 10,
    legend_width_mm = 14,
    
    top_margin_mm = 8,
    bottom_margin_mm = 8,
    left_margin_mm = 8,
    right_margin_mm = 8
) {
  
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(grid)
  
  # ==================== Dane ====================
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
      summary_rg   = first(na.omit(summary_rg)),
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
  ) %>%
    as_tibble()
  
  pairwise_rg_full <- full_grid %>%
    left_join(pairwise_sym, by = c("p1_id", "p2_id")) %>%
    mutate(
      summary_rg   = ifelse(p1_id == p2_id, 1, summary_rg),
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
  
  # ==================== Klastrowanie ====================
  rg_for_clust <- rg_matrix
  rg_for_clust[is.na(rg_for_clust)] <- 0
  
  if (cluster) {
    hc <- hclust(dist(rg_for_clust))
    trait_order <- hc$labels[hc$order]
  } else {
    trait_order <- selected_traits
  }
  
  rg_matrix <- rg_matrix[trait_order, trait_order, drop = FALSE]
  p_matrix  <- p_matrix[trait_order, trait_order, drop = FALSE]
  
  # ==================== Etykiety ====================
  label_vector <- metadata_ieu_EURsampleSize10000 %>%
    select(id, trait) %>%
    distinct() %>%
    mutate(label = paste0(id, " | ", trait)) %>%
    { setNames(.$label, .$id) }
  
  missing_ids <- setdiff(rownames(rg_matrix), names(label_vector))
  label_vector[missing_ids] <- missing_ids
  
  row_labels <- label_vector[rownames(rg_matrix)]
  col_labels <- label_vector[colnames(rg_matrix)]
  
  # ==================== Gwiazdki ====================
  get_stars <- function(p) {
    ifelse(
      is.na(p), "",
      ifelse(p < 0.001, "***",
             ifelse(p < 0.01, "**",
                    ifelse(p < 0.05, "*", "")))
    )
  }
  
  stars_matrix <- apply(p_matrix, c(1, 2), get_stars)
  
  # ==================== Rozmiary tekstu ====================
  max_row_w <- max_text_width(row_labels, gp = gpar(fontsize = label_size))
  max_col_h <- max_text_width(col_labels, gp = gpar(fontsize = label_size))
  # dla kolumn obróconych o 90 stopni szerokość tekstu staje się wysokością
  
  row_names_max_width    <- min(max_row_w + unit(3, "mm"), unit(95, "mm"))
  column_names_max_height <- min(max_col_h + unit(4, "mm"), unit(75, "mm"))
  
  # ==================== Kolory ====================
  col_fun <- colorRamp2(
    c(-1, 0, 1),
    c(low_color, mid_color, high_color)
  )
  
  # ==================== Tekst w komórkach ====================
  cell_fun <- function(j, i, x, y, width, height, fill) {
    if (!is.na(rg_matrix[i, j])) {
      grid.text(
        sprintf("%.2f", rg_matrix[i, j]),
        x = x,
        y = y + unit(1.2, "mm"),
        gp = gpar(fontsize = rg_size, fontface = "bold")
      )
    }
    
    if (stars_matrix[i, j] != "") {
      grid.text(
        stars_matrix[i, j],
        x = x,
        y = y - unit(2.0, "mm"),
        gp = gpar(fontsize = star_size, col = "#1a1a1a", fontface = "bold")
      )
    }
  }
  
  # ==================== Heatmapa ====================
  ht <- Heatmap(
    rg_matrix,
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
    
    show_heatmap_legend = FALSE,
    
    column_title = title,
    column_title_gp = gpar(fontsize = 12, fontface = "bold")
  )
  
  # ==================== Osobna legenda ====================
  lgd <- Legend(
    title = "Genetic\ncorrelation (rg)",
    col_fun = col_fun,
    at = c(-1, -0.5, 0, 0.5, 1),
    direction = "vertical",
    legend_height = unit(55, "mm"),
    grid_width = unit(6, "mm"),
    labels_gp = gpar(fontsize = 8),
    title_gp = gpar(fontsize = 9, fontface = "bold")
  )
  
  # ==================== Layout całej strony ====================
  total_right_block <- legend_gap_mm + legend_width_mm
  
  grid.newpage()
  pushViewport(
    viewport(
      layout = grid.layout(
        nrow = 1,
        ncol = 2,
        widths = unit.c(
          unit(1, "npc") - unit(total_right_block + right_margin_mm + left_margin_mm, "mm"),
          unit(total_right_block, "mm")
        )
      )
    )
  )
  
  # heatmapa
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  draw(
    ht,
    newpage = FALSE,
    padding = unit(c(top_margin_mm, 2, bottom_margin_mm, left_margin_mm), "mm")
  )
  popViewport()
  
  # legenda
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
  pushViewport(
    viewport(
      x = unit(legend_gap_mm / total_right_block, "npc"),
      y = unit(0.5, "npc"),
      just = c("left", "center"),
      width = unit(legend_width_mm, "mm"),
      height = unit(1, "npc") - unit(top_margin_mm + bottom_margin_mm, "mm")
    )
  )
  draw(lgd, just = c("left", "center"))
  popViewport(2)
  
  invisible(list(heatmap = ht, legend = lgd))
}

create_rg_heatmap(
  selected_traits = selected_traits,
  ldsc_results_annotated = ldsc_results_annotated,
  metadata_ieu_EURsampleSize10000 = metadata_ieu_EURsampleSize10000,
  label_size = 8.5,
  rg_size = 7.5,
  star_size = 9,
  legend_gap_mm = -10
)
