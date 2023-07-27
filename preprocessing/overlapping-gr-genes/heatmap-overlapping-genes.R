
gr_gene_database_preprocessing %>% 
  filter(trimws(Gene_name) != "") -> gr_gene_database_preprocessing

gr_gene_database_preprocessing %>% 
  filter(label == "34362910_NA_NA")


pairwise_gene_overlap(gr_gene_database_preprocessing, "Gene_name", "label") -> overlap_matrix

gene_set_sizes <- table(unlist(gr_gene_database_preprocessing$label))

normalized_overlap_matrix <- sweep(overlap_matrix, MARGIN = 1, STATS = gene_set_sizes, FUN = "/")

color_palette <- colorRampPalette(c("white", "red"))(100)


pheatmap(overlap_matrix, color = color_palette, cluster_rows = TRUE, cluster_cols = TRUE)


pheatmap(normalized_overlap_matrix, color = color_palette, cluster_rows = FALSE, cluster_cols = FALSE)

Heatmap(normalized_ordered, cluster_rows = F, cluster_columns = F)

# Create the heatmap
Heatmap(overlap_matrix)

Heatmap(overlap_matrix, 
        name = "overlap", 
        col = color_palette, 
        show_row_names = TRUE, 
        show_column_names = TRUE)


test_result <- chisq.test(overlap_matrix)

# print the result
print(test_result)

# perform the Fisher's Exact Test
test_result <- fisher.test(overlap_matrix, simulate.p.value = TRUE)

# print the result
print(test_result)
