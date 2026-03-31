# ================================================
# FINALNA FUNKCJA – Genetic Correlation Heatmap
# ================================================

create_rg_heatmap <- function(
    selected_traits,
    ldsc_results_annotated,
    metadata_ieu_EURsampleSize10000,
    
    # Tytuł i podstawowe opcje
    title = "Genetic correlation heatmap (LDSC)",
    
    # Rozmiary czcionek
    label_size = 10,      # etykiety traitów
    rg_size = 9.2,        # wartości rg
    star_size = 12,       # gwiazdki
    
    # Kolory gradientu
    low_color  = "#2166AC",
    mid_color  = "#F7F7F7",
    high_color = "#B2182B",
    
    # Położenie etykiet
    row_labels_side = "left",   # "left" lub "right"
    col_labels_side = "top",    # "top" lub "bottom"
    
    # Klastrowanie
    cluster = TRUE
) {
  
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(tidyr)
  
  # ==================== Przygotowanie danych ====================
  
  pairwise_rg <- ldsc_results_annotated %>%
    filter(p1_id %in% selected_traits, p2_id %in% selected_traits) %>%
    filter(!is.na(summary_rg)) %>%
    select(p1_id, p2_id, summary_rg, summary_rg_p)
  
  # Symetryzacja par
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
    pairwise_unique %>% transmute(p1_id = trait_min, p2_id = trait_max, summary_rg, summary_rg_p),
    pairwise_unique %>% transmute(p1_id = trait_max, p2_id = trait_min, summary_rg, summary_rg_p)
  )
  
  # Pełna siatka
  full_grid <- expand.grid(p1_id = selected_traits, 
                           p2_id = selected_traits, 
                           stringsAsFactors = FALSE) %>%
    as_tibble()
  
  pairwise_rg_full <- full_grid %>%
    left_join(pairwise_sym, by = c("p1_id", "p2_id")) %>%
    mutate(
      summary_rg = ifelse(p1_id == p2_id, 1, summary_rg),
      summary_rg_p = ifelse(p1_id == p2_id, NA, summary_rg_p)
    )
  
  # Macierze
  rg_matrix <- pairwise_rg_full %>%
    select(p1_id, p2_id, summary_rg) %>%
    pivot_wider(names_from = p2_id, values_from = summary_rg) %>%
    column_to_rownames("p1_id") %>%
    as.matrix()
  
  p_matrix <- pairwise_rg_full %>%
    select(p1_id, p2_id, summary_rg_p) %>%
    pivot_wider(names_from = p2_id, values_from = summary_rg_p) %>%
    column_to_rownames("p1_id") %>%
    as.matrix()
  
  # Klastrowanie
  rg_matrix_for_clust <- rg_matrix
  rg_matrix_for_clust[is.na(rg_matrix_for_clust)] <- 0
  
  hc <- hclust(dist(rg_matrix_for_clust))
  trait_order <- hc$labels[hc$order]
  
  rg_matrix <- rg_matrix[trait_order, trait_order]
  p_matrix  <- p_matrix[trait_order, trait_order]
  
  # Etykiety
  trait_labels_df <- metadata_ieu_EURsampleSize10000 %>%
    select(id, trait) %>%
    distinct() %>%
    mutate(label = paste0(id, " | ", trait))
  
  label_vector <- trait_labels_df$label
  names(label_vector) <- trait_labels_df$id
  
  missing_ids <- setdiff(rownames(rg_matrix), names(label_vector))
  label_vector[missing_ids] <- missing_ids
  
  # Gwiazdkki
  get_stars <- function(p) {
    ifelse(is.na(p), "",
           ifelse(p < 0.001, "***",
                  ifelse(p < 0.01, "**",
                         ifelse(p < 0.05, "*", ""))))
  }
  
  stars_matrix <- apply(p_matrix, c(1,2), get_stars)
  
  # Kolory
  col_fun <- colorRamp2(c(-1, 0, 1), c(low_color, mid_color, high_color))
  
  # cell_fun – Twoje ustawienia
  cell_fun <- function(j, i, x, y, width, height, fill) {
    rg_val <- rg_matrix[i, j]
    
    if (!is.na(rg_val)) {
      grid.text(sprintf("%.2f", rg_val),
                x, y - unit(-1, "mm"),
                gp = gpar(fontsize = rg_size, fontface = "bold"))
    }
    
    star <- stars_matrix[i, j]
    if (star != "") {
      grid.text(star,
                x, y + unit(-3, "mm"),
                gp = gpar(fontsize = star_size, col = "#1a1a1a", fontface = "bold"))
    }
  }
  
  # ==================== Heatmapa ====================
  
  ht <- Heatmap(
    rg_matrix,
    name = "rg",
    
    col = col_fun,
    rect_gp = gpar(col = "white", lwd = 1),
    
    cluster_rows = cluster,
    cluster_columns = cluster,
    show_row_dend = TRUE,
    show_column_dend = TRUE,
    row_dend_width = unit(3.8, "cm"),
    column_dend_height = unit(3.8, "cm"),
    
    # <<< NOWE: Położenie etykiet >>>
    row_labels = label_vector[rownames(rg_matrix)],
    column_labels = label_vector[colnames(rg_matrix)],
    
    row_names_side = row_labels_side,     # "left" lub "right"
    column_names_side = col_labels_side,  # "top" lub "bottom"
    
    row_names_gp = gpar(fontsize = label_size),
    column_names_gp = gpar(fontsize = label_size),
    column_names_rot = 90,
    
    cell_fun = cell_fun,
    
    heatmap_legend_param = list(
      title = "Genetic\ncorrelation (rg)",
      at = c(-1, -0.5, 0, 0.5, 1)
    ),
    
    column_title = title,
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    
    width = unit(0.92, "npc"),
    height = unit(0.92, "npc")
  )
  
  draw(ht, padding = unit(c(12, 40, 12, 40), "mm"))
  
  invisible(ht)
}


create_rg_heatmap(
  selected_traits = selected_traits,
  ldsc_results_annotated = ldsc_results_annotated,
  metadata_ieu_EURsampleSize10000 = metadata_ieu_EURsampleSize10000,
  
  row_labels_side = "left",   # ← nowe
  col_labels_side = "top",  # ← nowe
  
  label_size = 9.5,
  rg_size = 9,
  star_size = 13
)
