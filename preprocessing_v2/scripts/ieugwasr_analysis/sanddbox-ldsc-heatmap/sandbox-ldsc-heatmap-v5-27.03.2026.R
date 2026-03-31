library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

# =========================================================
# 1. wybór traitów powiązanych z traitem referencyjnym
# =========================================================
reference_trait <- "ieu-a-806"

selected_traits <- ldsc_results_annotated %>%
  filter(p1_id == reference_trait, summary_rg_p < 0.05) %>%
  pull(p2_id) %>%
  unique()

selected_traits <- unique(c(reference_trait, selected_traits))

# =========================================================
# 2. pobranie wszystkich wyników dla wybranych traitów
# =========================================================
pairwise_rg <- ldsc_results_annotated %>%
  filter(p1_id %in% selected_traits, p2_id %in% selected_traits) %>%
  select(p1_id, p2_id, summary_rg, summary_rg_p)

# =========================================================
# 3. sprowadzenie A-B i B-A do jednej wspólnej pary
# =========================================================
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

# =========================================================
# 4. rozpisanie z powrotem na obie strony
# =========================================================
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

# =========================================================
# 5. pełna siatka wszystkich kombinacji
# =========================================================
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

# =========================================================
# 6. etykiety: id | trait
# =========================================================
trait_labels_df <- metadata_ieu_EURsampleSize10000 %>%
  select(id, trait) %>%
  distinct() %>%
  mutate(label = paste0(id, " | ", trait))

label_vector <- trait_labels_df$label
names(label_vector) <- trait_labels_df$id

missing_ids <- setdiff(selected_traits, names(label_vector))
label_vector[missing_ids] <- missing_ids

# =========================================================
# 7. macierz rg
# =========================================================
rg_matrix <- pairwise_rg_full %>%
  select(p1_id, p2_id, summary_rg) %>%
  pivot_wider(names_from = p2_id, values_from = summary_rg) %>%
  as.data.frame()

rownames(rg_matrix) <- rg_matrix$p1_id
rg_matrix$p1_id <- NULL
rg_matrix <- as.matrix(rg_matrix)

# =========================================================
# 8. macierz p-value
# =========================================================
p_matrix <- pairwise_rg_full %>%
  select(p1_id, p2_id, summary_rg_p) %>%
  pivot_wider(names_from = p2_id, values_from = summary_rg_p) %>%
  as.data.frame()

rownames(p_matrix) <- p_matrix$p1_id
p_matrix$p1_id <- NULL
p_matrix <- as.matrix(p_matrix)

# =========================================================
# 9. upewnienie się co do kolejności
# =========================================================
rg_matrix <- rg_matrix[selected_traits, selected_traits]
p_matrix  <- p_matrix[selected_traits, selected_traits]

# =========================================================
# 10. klastrowanie bez absolute values
#     rg traktujemy jako similarity
#     distance = (1 - rg)/2
# =========================================================
rg_matrix_for_clustering <- rg_matrix
rg_matrix_for_clustering[is.na(rg_matrix_for_clustering)] <- 0

distance_matrix <- (1 - rg_matrix_for_clustering) / 2

# bezpieczeństwo numeryczne
distance_matrix <- (distance_matrix + t(distance_matrix)) / 2
diag(distance_matrix) <- 0

hc <- hclust(as.dist(distance_matrix), method = "average")

ordered_traits <- hc$labels[hc$order]

rg_matrix <- rg_matrix[ordered_traits, ordered_traits]
p_matrix  <- p_matrix[ordered_traits, ordered_traits]

# =========================================================
# 11. etykiety po uporządkowaniu
# =========================================================
row_labels <- label_vector[rownames(rg_matrix)]
col_labels <- label_vector[colnames(rg_matrix)]

# =========================================================
# 12. funkcje pomocnicze
# =========================================================
p_to_stars <- function(p) {
  if (is.na(p)) {
    return("")
  } else if (p < 0.001) {
    return("***")
  } else if (p < 0.01) {
    return("**")
  } else if (p < 0.05) {
    return("*")
  } else {
    return("")
  }
}

format_rg <- function(x) {
  if (is.na(x)) {
    return("")
  } else {
    return(sprintf("%.2f", x))
  }
}

# macierz gwiazdek
star_matrix <- matrix(
  vapply(as.vector(p_matrix), p_to_stars, character(1)),
  nrow = nrow(p_matrix),
  ncol = ncol(p_matrix),
  byrow = FALSE
)

# macierz tekstowa z wartościami rg
rg_label_matrix <- matrix(
  vapply(as.vector(rg_matrix), format_rg, character(1)),
  nrow = nrow(rg_matrix),
  ncol = ncol(rg_matrix),
  byrow = FALSE
)

# =========================================================
# 13. skala kolorów
# =========================================================
col_fun <- colorRamp2(
  c(-1, 0, 1),
  c("blue", "white", "red")
)

# =========================================================
# 14. heatmapa
# =========================================================
ht <- Heatmap(
  rg_matrix,
  name = "rg",
  col = col_fun,
  na_col = "grey85",
  
  cluster_rows = as.dendrogram(hc),
  cluster_columns = as.dendrogram(hc),
  
  row_labels = row_labels,
  column_labels = col_labels,
  
  row_names_side = "left",
  column_names_side = "bottom",
  column_names_rot = 90,
  
  rect_gp = gpar(col = "white", lwd = 1),
  
  heatmap_legend_param = list(
    title = "rg",
    at = c(-1, -0.5, 0, 0.5, 1)
  ),
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    rg_lab <- rg_label_matrix[i, j]
    star_lab <- star_matrix[i, j]
    
    if (rg_lab != "") {
      grid.text(
        rg_lab,
        x = x,
        y = y + unit(2, "pt"),
        gp = gpar(fontsize = 6)
      )
    }
    
    if (star_lab != "") {
      grid.text(
        star_lab,
        x = x,
        y = y - unit(6, "pt"),
        gp = gpar(fontsize = 6, fontface = "bold")
      )
    }
  },
  
  column_title = "Genetic correlation heatmap",
  row_title = NULL
)

draw(ht)
