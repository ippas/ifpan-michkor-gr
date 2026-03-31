library(ComplexHeatmap)
library(circlize)

# =========================================
# Macierze (zakładam, że masz już rg_matrix i p_matrix w kolejności trait_order)
# =========================================
rg_matrix <- rg_matrix[trait_order, trait_order]
p_matrix  <- p_matrix[trait_order, trait_order]

# Gwiazdkki istotności
get_stars <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01,  "**",
                       ifelse(p < 0.05,  "*", ""))))
}

stars_matrix <- apply(p_matrix, c(1,2), get_stars)

# Kolory
col_fun <- colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B"))  # ładniejszy niebiesko-czerwony

# =========================================
# cell_fun – rg wyżej, gwiazdki niżej
# =========================================
cell_fun <- function(j, i, x, y, width, height, fill) {
  rg_val <- rg_matrix[i, j]
  
  if (!is.na(rg_val)) {
    # rg wyżej i pogrubiona
    grid.text(sprintf("%.2f", rg_val), 
              x, y - unit(3.5, "mm"),
              gp = gpar(fontsize = 9.5, fontface = "bold"))
  }
  
  # gwiazdki niżej
  star <- stars_matrix[i, j]
  if (star != "") {
    grid.text(star, 
              x, y + unit(4, "mm"),
              gp = gpar(fontsize = 12, col = "black", fontface = "bold"))
  }
}

# =========================================
# Heatmapa z dendrogramami
# =========================================
ht <- Heatmap(
  rg_matrix,
  name = "rg",
  
  col = col_fun,
  rect_gp = gpar(col = "white", lwd = 0.8),
  
  # === Dendrogramy ===
  cluster_rows = TRUE,           # włączamy klastrowanie (użyje hc z Twojego kodu)
  cluster_columns = TRUE,
  show_row_dend = TRUE,          # dendrogram po lewej
  show_column_dend = TRUE,       # dendrogram na górze
  row_dend_width = unit(4, "cm"), 
  column_dend_height = unit(4, "cm"),
  
  # === Etykiety ===
  row_labels = label_vector[rownames(rg_matrix)],
  column_labels = label_vector[colnames(rg_matrix)],
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  column_names_rot = 90,         # podpisy na górze pionowo
  
  # cell_fun
  cell_fun = cell_fun,
  
  # Legenda
  heatmap_legend_param = list(
    title = "Genetic correlation (rg)",
    at = c(-1, -0.5, 0, 0.5, 1),
    legend_height = unit(6, "cm")
  ),
  
  column_title = "Genetic correlation heatmap (LDSC)",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  
  width = unit(0.85, "npc"),
  height = unit(0.85, "npc")
)

# Rysowanie
draw(ht, 
     padding = unit(c(10, 25, 10, 25), "mm"),   # miejsce na dendrogramy i etykiety
     merge_legend = TRUE)