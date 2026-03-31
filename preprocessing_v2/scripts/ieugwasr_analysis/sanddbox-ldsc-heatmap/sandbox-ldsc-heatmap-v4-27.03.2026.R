library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

# =========================================
# 1. wybór traitów powiązanych z ieu-a-806
# =========================================
selected_traits <- ldsc_results_annotated %>%
  filter(p1_id == "ieu-a-806", summary_rg_p < 0.05) %>%
  pull(p2_id) %>%
  unique()

selected_traits <- unique(c("ieu-a-806", selected_traits))

# =========================================
# 2. pobranie wszystkich wyników dla wybranych traitów
# =========================================
pairwise_rg <- ldsc_results_annotated %>%
  filter(p1_id %in% selected_traits, p2_id %in% selected_traits) %>%
  filter(!is.na(summary_rg)) %>%
  select(p1_id, p2_id, summary_rg, summary_rg_p)

# =========================================
# 3. sprowadzenie A-B i B-A do jednej wspólnej pary
# =========================================
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

# =========================================
# 4. rozpisanie z powrotem na obie strony
# =========================================
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

# =========================================
# 5. pełna siatka wszystkich kombinacji
# =========================================
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
# 6. etykiety: id | trait
# =========================================
trait_labels_df <- metadata_ieu_EURsampleSize10000 %>%
  select(id, trait) %>%
  distinct() %>%
  mutate(label = paste0(id, " | ", trait))

label_vector <- trait_labels_df$label
names(label_vector) <- trait_labels_df$id

missing_ids <- setdiff(selected_traits, names(label_vector))
label_vector[missing_ids] <- missing_ids

# =========================================
# 7. macierze rg i p
# =========================================
rg_matrix <- pairwise_rg_full %>%
  select(p1_id, p2_id, summary_rg) %>%
  pivot_wider(names_from = p2_id, values_from = summary_rg) %>%
  as.data.frame()

rownames(rg_matrix) <- rg_matrix$p1_id
rg_matrix$p1_id <- NULL
rg_matrix <- as.matrix(rg_matrix)

p_matrix <- pairwise_rg_full %>%
  select(p1_id, p2_id, summary_rg_p) %>%
  pivot_wider(names_from = p2_id, values_from = summary_rg_p) %>%
  as.data.frame()

rownames(p_matrix) <- p_matrix$p1_id
p_matrix$p1_id <- NULL
p_matrix <- as.matrix(p_matrix)

rg_matrix <- rg_matrix[selected_traits, selected_traits]
p_matrix  <- p_matrix[selected_traits, selected_traits]

# =========================================
# 8. klastrowanie
#    chcemy grupować dodatnie z dodatnimi
#    i ujemne z ujemnymi
# =========================================
rg_matrix_for_clustering <- rg_matrix
rg_matrix_for_clustering[is.na(rg_matrix_for_clustering)] <- 0

# WARIANT 1: prosty i zwykle bardzo dobry
hc <- hclust(dist(rg_matrix_for_clustering), method = "complete")

# alternatywnie możesz testować:
# hc <- hclust(dist(rg_matrix_for_clustering), method = "ward.D2")

ordered_traits <- hc$labels[hc$order]

rg_matrix <- rg_matrix[ordered_traits, ordered_traits]
p_matrix  <- p_matrix[ordered_traits, ordered_traits]

# =========================================
# 9. nazwy osi
# =========================================
row_labels <- label_vector[rownames(rg_matrix)]
col_labels <- label_vector[colnames(rg_matrix)]

# =========================================
# 10. kolory
# =========================================
col_fun <- colorRamp2(
  c(-1, 0, 1),
  c("blue", "white", "red")
)

# =========================================
# 11. funkcja do gwiazdek
# =========================================
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

star_matrix <- matrix(
  vapply(as.vector(p_matrix), p_to_stars, character(1)),
  nrow = nrow(p_matrix),
  ncol = ncol(p_matrix),
  byrow = FALSE
)

# =========================================
# 12. heatmapa
# =========================================
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
    lab <- star_matrix[i, j]
    if (lab != "") {
      grid.text(lab, x, y, gp = gpar(fontsize = 7))
    }
  },
  
  column_title = "Genetic correlation heatmap",
  row_title = NULL
)

draw(ht)