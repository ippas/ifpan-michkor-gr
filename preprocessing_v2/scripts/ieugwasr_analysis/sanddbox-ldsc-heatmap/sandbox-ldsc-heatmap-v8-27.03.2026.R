library(ComplexHeatmap)
library(circlize)

# Kolory (ładniejszy kontrast)
col_fun <- colorRamp2(c(-1, 0, 1), c("#2166AC", "#F7F7F7", "#B2182B"))

# cell_fun z lepszym wyśrodkowaniem
cell_fun <- function(j, i, x, y, width, height, fill) {
  rg_val <- rg_matrix[i, j]
  
  if (!is.na(rg_val)) {
    grid.text(sprintf("%.2f", rg_val), 
              x, y - unit(-1, "mm"),           # bliżej środka
              gp = gpar(fontsize = 9.2, fontface = "bold"))
  }
  
  star <- stars_matrix[i, j]
  if (star != "") {
    grid.text(star, 
              x, y + unit(-3, "mm"),           # gwiazdki niżej
              gp = gpar(fontsize = 12, col = "#1a1a1a", fontface = "bold"))
  }
}

# Heatmapa
ht <- Heatmap(
  rg_matrix,
  name = "rg",
  
  col = col_fun,
  rect_gp = gpar(col = "white", lwd = 1),
  
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_dend = TRUE,
  show_column_dend = TRUE,
  row_dend_width = unit(3.5, "cm"),
  column_dend_height = unit(3.5, "cm"),
  
  row_labels = label_vector[rownames(rg_matrix)],
  column_labels = label_vector[colnames(rg_matrix)],
  
  row_names_gp = gpar(fontsize = 9.5),
  column_names_gp = gpar(fontsize = 9.5),
  column_names_rot = 90,
  
  cell_fun = cell_fun,
  
  heatmap_legend_param = list(
    title = "Genetic\ncorrelation (rg)",
    at = c(-1, -0.5, 0, 0.5, 1)
  ),
  
  column_title = "Genetic correlation heatmap (LDSC)",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  
  width = unit(0.9, "npc"),
  height = unit(0.9, "npc")
)

draw(ht, padding = unit(c(8, 30, 8, 30), "mm"))

