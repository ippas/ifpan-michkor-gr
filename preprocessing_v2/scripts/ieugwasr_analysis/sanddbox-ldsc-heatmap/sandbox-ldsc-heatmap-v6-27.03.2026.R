library(ComplexHeatmap)
library(circlize)   # do colorRamp2

# =========================================
# Przygotowanie macierzy wartości rg i p-wartości
# =========================================

# Macierz z rg (już masz)
rg_matrix <- pairwise_rg_full %>%
  select(p1_id, p2_id, summary_rg) %>%
  pivot_wider(names_from = p2_id, values_from = summary_rg) %>%
  column_to_rownames("p1_id") %>%
  as.matrix()

# Macierz z p-wartościami (analogicznie)
p_matrix <- pairwise_rg_full %>%
  select(p1_id, p2_id, summary_rg_p) %>%
  pivot_wider(names_from = p2_id, values_from = summary_rg_p) %>%
  column_to_rownames("p1_id") %>%
  as.matrix()

# Kolejność z Twojego klastrowania (zachowujemy to samo grupowanie)
rg_matrix <- rg_matrix[trait_order, trait_order]
p_matrix  <- p_matrix[trait_order, trait_order]

# =========================================
# Funkcja do gwiazdek istotności
# =========================================
get_stars <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01,  "**",
                       ifelse(p < 0.05,  "*", ""))))
}

stars_matrix <- apply(p_matrix, c(1,2), get_stars)

# =========================================
# Kolory dla rg (ładniejszy gradient niż w ggplot)
# =========================================
col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

# =========================================
# cell_fun – najważniejsze: wartości rg + gwiazdki pod spodem
# =========================================
cell_fun <- function(j, i, x, y, width, height, fill) {
  rg_val <- rg_matrix[i, j]
  
  # Wartość rg (zaokrąglona do 2 miejsc)
  if (!is.na(rg_val)) {
    grid.text(sprintf("%.2f", rg_val), 
              x, y - unit(2, "mm"),   # trochę wyżej
              gp = gpar(fontsize = 9, fontface = "bold"))
  }
  
  # Gwiazdki istotności (niżej)
  star <- stars_matrix[i, j]
  if (star != "") {
    grid.text(star, 
              x, y + unit(3, "mm"),   # trochę niżej
              gp = gpar(fontsize = 11, col = "black", fontface = "bold"))
  }
}

# =========================================
# Rysowanie heatmapy
# =========================================
ht <- Heatmap(
  rg_matrix,
  name = "rg",
  
  col = col_fun,
  rect_gp = gpar(col = "white", lwd = 1),   # białe linie między kafelkami
  
  # Etykiety z Twojego label_vector
  row_labels = label_vector[rownames(rg_matrix)],
  column_labels = label_vector[colnames(rg_matrix)],
  
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  column_names_rot = 90,
  
  cluster_rows = FALSE,     # już masz własne klastrowanie
  cluster_columns = FALSE,
  
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(title = "Genetic\ncorrelation (rg)",
                              at = c(-1, -0.5, 0, 0.5, 1)),
  
  cell_fun = cell_fun,
  
  # Opcjonalnie: tytuł i marginesy
  column_title = "Genetic correlation heatmap (LDSC)",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  
  width = unit(1, "npc") * 0.9,   # żeby nie był za szeroki
  height = unit(1, "npc") * 0.9
)

# Rysuj
draw(ht, padding = unit(c(2, 20, 2, 20), "mm"))   # trochę miejsca na etykiety