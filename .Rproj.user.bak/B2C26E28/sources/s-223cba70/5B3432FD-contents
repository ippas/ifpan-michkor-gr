# Define the custom color palette
custom_palette <- c(
  "pastel_blue"       = "white",
  "pastel_light_blue" = "#f8dedd",
  "white"             = "#f1bcbb",
  "pastel_orange"     = "#edacab",
  "pastel_red"        = "#e68a89"
)

# Function to create and plot a triangular heatmap
plot_triangular_heatmap <- function(mat, palette) {
  
  # Compute the Euclidean distance matrix for rows
  row_dist <- dist(mat, method = "euclidean")
  
  # Perform hierarchical clustering on rows using the "complete" linkage method
  row_hclust <- hclust(row_dist, method = "complete")
  
  # Compute the Euclidean distance matrix for columns
  col_dist <- dist(t(mat), method = "euclidean")
  
  # Perform hierarchical clustering on columns using the "complete" linkage method
  col_hclust <- hclust(col_dist, method = "complete")
  
  # Reorder matrix based on clustering
  reordered_mat <- mat[row_hclust$order, col_hclust$order]
  
  reordered_mat <- mat
  reordered_mat <-  mat[nrow(mat):1, ncol(mat):1]

  
  # Create and plot the heatmap using ComplexHeatmap with custom cell_fun
  ht <- Heatmap(
    reordered_mat,
    name = "Jaccard Index",
    na_col = "white",
    col = palette,
    show_row_names = TRUE,
    show_column_names = TRUE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    rect_gp = gpar(type = "none"),
    row_names_side = "left",
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (i >= j) {
        grid.rect(x = x, y = y, width = width, height = height,
                  gp = gpar(fill = fill, col = fill))
        grid.text(sprintf("%.2f", reordered_mat[i, j]), x = x, y = y, 
                  gp = gpar(fontsize = 8, col = "black"))
      }
    },
    heatmap_legend_param = list(direction = "horizontal")
  )
  draw(ht, heatmap_legend_side = "bottom")
}


tissue_master_df %>% 
  filter(regulation == "down", 
         rank_criterion == "log2ratio",
         size_list == 50) %>% 
  ungroup(hgnc_symbol) %>%
  select(c(tissue, hgnc_symbol)) %>% 
  as.data.frame() %>% 
  split(.$tissue, .$hgnc_symbol ) %>% lapply(., function(x) {x$hgnc_symbol}) %>%  
  append(list(master_up_double_rs = {regulation_master_double_rank_score_lists %>% 
      bind_rows(.id = "label") %>% 
      separate(label, into = c("regulation", "metric", "size_list")) %>% 
      filter(metric == "log2ratio") %>% 
      filter(regulation == "down") %>% 
      .$hgnc_symbol})) %>% 
  jaccard_matrix(gene_list = .) %>% plot_triangular_heatmap(mat = ., as.vector(custom_palette))

tissue_master_df %>% 
  filter(regulation == "up", 
         rank_criterion == "log2ratio",
         size_list == 50) %>% 
  ungroup(hgnc_symbol) %>%
  select(c(tissue, hgnc_symbol)) %>% 
  as.data.frame() %>% 
  split(.$tissue, .$hgnc_symbol ) %>% lapply(., function(x) {x$hgnc_symbol}) %>%  
  append(list(master_up_double_rs = {regulation_master_double_rank_score_lists %>% 
      bind_rows(.id = "label") %>% 
      separate(label, into = c("regulation", "metric", "size_list")) %>% 
      filter(metric == "log2ratio") %>% 
      filter(regulation == "up") %>% 
      .$hgnc_symbol})) %>% 
  jaccard_matrix(gene_list = .) %>% plot_triangular_heatmap(mat = ., as.vector(custom_palette))




################################################################################
tissue_master_df %>% 
  filter(regulation == "down", 
         rank_criterion == "log2ratio",
         size_list == 50) %>% 
  ungroup(hgnc_symbol) %>%
  select(c(tissue, hgnc_symbol)) %>% 
  as.data.frame() %>% 
  split(.$tissue, .$hgnc_symbol ) %>% lapply(., function(x) {x$hgnc_symbol}) %>%  
  append(list(master_up_double_rs = {regulation_master_double_rank_score_lists %>% 
      bind_rows(.id = "label") %>% 
      separate(label, into = c("regulation", "metric", "size_list")) %>% 
      filter(metric == "log2ratio") %>% 
      filter(regulation == "down") %>% 
      .$hgnc_symbol})) %>% gene_list_to_matrix() %>% 
  t() %>%
  proxy::dist(method = "Jaccard") %>%
  as.matrix() %>% { 1 - . }  %>% pheatmap(display_numbers = T)
