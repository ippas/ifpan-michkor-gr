# ================================================
# WERSJA Z LEPSZYM ODSTĘPEM OD LEGENDY + MNIEJ BIAŁEJ PRZESTRZENI
# ================================================

create_rg_heatmap <- function(
    selected_traits,
    ldsc_results_annotated,
    metadata_ieu_EURsampleSize10000,
    
    title = "Genetic correlation heatmap (LDSC)",
    
    label_size = 9.8,
    rg_size = 9.2,
    star_size = 12,
    
    low_color  = "#2166AC",
    mid_color  = "#F7F7F7",
    high_color = "#B2182B",
    
    row_labels_side = "left",
    col_labels_side = "top",
    
    cluster = TRUE,
    auto_adjust = TRUE,
    
    # NOWY parametr – kontroluje odstęp między heatmapą a legendą
    legend_gap_mm = 18      # zwiększ jeśli chcesz większą przerwę (np. 25-30)
) {
  
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(tidyr)
  library(grid)
  
  # ==================== Przygotowanie macierzy ====================
  # (ta sama sekcja co wcześniej – wklejam skróconą, ale pełną wersję musisz mieć całą)
  
  pairwise_rg <- ldsc_results_annotated %>%
    filter(p1_id %in% selected_traits, p2_id %in% selected_traits) %>%
    filter(!is.na(summary_rg)) %>%
    select(p1_id, p2_id, summary_rg, summary_rg_p)
  
  pairwise_unique <- pairwise_rg %>%
    rowwise() %>%
    mutate(trait_min = sort(c(p1_id, p2_id))[1],
           trait_max = sort(c(p1_id, p2_id))[2]) %>%
    ungroup() %>%
    group_by(trait_min, trait_max) %>%
    summarise(summary_rg = first(na.omit(summary_rg)),
              summary_rg_p = first(na.omit(summary_rg_p)), .groups = "drop")
  
  pairwise_sym <- bind_rows(
    pairwise_unique %>% transmute(p1_id = trait_min, p2_id = trait_max, summary_rg, summary_rg_p),
    pairwise_unique %>% transmute(p1_id = trait_max, p2_id = trait_min, summary_rg, summary_rg_p)
  )
  
  full_grid <- expand.grid(p1_id = selected_traits, p2_id = selected_traits, stringsAsFactors = FALSE) %>% as_tibble()
  
  pairwise_rg_full <- full_grid %>%
    left_join(pairwise_sym, by = c("p1_id", "p2_id")) %>%
    mutate(summary_rg = ifelse(p1_id == p2_id, 1, summary_rg),
           summary_rg_p = ifelse(p1_id == p2_id, NA, summary_rg_p))
  
  rg_matrix <- pairwise_rg_full %>%
    pivot_wider(id_cols = p1_id, names_from = p2_id, values_from = summary_rg) %>%
    column_to_rownames("p1_id") %>% as.matrix()
  
  p_matrix <- pairwise_rg_full %>%
    pivot_wider(id_cols = p1_id, names_from = p2_id, values_from = summary_rg_p) %>%
    column_to_rownames("p1_id") %>% as.matrix()
  
  # Klastrowanie + kolejność
  rg_for_clust <- rg_matrix; rg_for_clust[is.na(rg_for_clust)] <- 0
  hc <- hclust(dist(rg_for_clust))
  trait_order <- hc$labels[hc$order]
  
  rg_matrix <- rg_matrix[trait_order, trait_order]
  p_matrix  <- p_matrix[trait_order, trait_order]
  
  # Etykiety
  label_vector <- metadata_ieu_EURsampleSize10000 %>%
    select(id, trait) %>% distinct() %>%
    mutate(label = paste0(id, " | ", trait)) %>%
    { setNames(.$label, .$id) }
  
  missing <- setdiff(rownames(rg_matrix), names(label_vector))
  label_vector[missing] <- missing
  
  # Gwiazdkki
  get_stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))))
  stars_matrix <- apply(p_matrix, c(1,2), get_stars)
  
  col_fun <- colorRamp2(c(-1, 0, 1), c(low_color, mid_color, high_color))
  
  cell_fun <- function(j, i, x, y, width, height, fill) {
    if (!is.na(rg_matrix[i,j])) {
      grid.text(sprintf("%.2f", rg_matrix[i,j]), x, y - unit(-1, "mm"),
                gp = gpar(fontsize = rg_size, fontface = "bold"))
    }
    if (stars_matrix[i,j] != "") {
      grid.text(stars_matrix[i,j], x, y + unit(-3, "mm"),
                gp = gpar(fontsize = star_size, col = "#1a1a1a", fontface = "bold"))
    }
  }
  
  # ==================== Automatyczne + ręczne odstępy ====================
  if (auto_adjust) {
    row_labels_text <- label_vector[rownames(rg_matrix)]
    max_row_w <- max_text_width(row_labels_text, gp = gpar(fontsize = label_size))
    
    dend_w_cm <- max(4.2, as.numeric(convertWidth(max_row_w + unit(10,"mm"), "cm")))
    
    padding_mm <- c(25, 95 + legend_gap_mm, 40, 25)   # right = 95 + legend_gap_mm
  } else {
    dend_w_cm <- 4.8
    padding_mm <- c(25, 110, 40, 25)
  }
  
  # Heatmapa
  ht <- Heatmap(
    rg_matrix,
    name = "rg",
    col = col_fun,
    rect_gp = gpar(col = "white", lwd = 1),
    
    cluster_rows = cluster,
    cluster_columns = cluster,
    show_row_dend = TRUE,
    show_column_dend = TRUE,
    
    row_dend_width = unit(dend_w_cm, "cm"),
    column_dend_height = unit(5.2, "cm"),
    
    row_names_max_width = unit(180, "mm"),   # ograniczenie, żeby nie było za szeroko
    
    row_labels = label_vector[rownames(rg_matrix)],
    column_labels = label_vector[colnames(rg_matrix)],
    
    row_names_side = row_labels_side,
    column_names_side = col_labels_side,
    
    row_names_gp = gpar(fontsize = label_size),
    column_names_gp = gpar(fontsize = label_size),
    column_names_rot = 90,
    
    cell_fun = cell_fun,
    
    heatmap_legend_param = list(
      title = "Genetic\ncorrelation (rg)",
      at = c(-1, -0.5, 0, 0.5, 1),
      legend_height = unit(7, "cm"),
      legend_width = unit(1.2, "cm")
    ),
    
    column_title = title,
    column_title_gp = gpar(fontsize = 14, fontface = "bold")
  )
  
  draw(ht, padding = unit(padding_mm, "mm"))
  
  invisible(ht)
}
create_rg_heatmap(
  selected_traits = selected_traits,
  ldsc_results_annotated = ldsc_results_annotated,
  metadata_ieu_EURsampleSize10000 = metadata_ieu_EURsampleSize10000,
  
  row_labels_side = "left",
  col_labels_side = "top",
  label_size = 9.5,
  legend_gap_mm = 60,
  auto_adjust = TRUE          # ← to jest klucz
)
