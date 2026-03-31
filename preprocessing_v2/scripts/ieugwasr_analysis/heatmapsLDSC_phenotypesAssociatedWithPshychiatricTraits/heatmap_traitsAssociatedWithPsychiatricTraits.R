metadata_ieu_EURsampleSize10000 %>% 
  filter(subcategory == "Psychiatric / neurological")


metadata_ieu_EURsampleSize10000 %>% 
  filter(subcategory == "Psychiatric / neurological") %>% 
  .$id -> psychiatric_id


selected_traits <- ldsc_results_annotated %>%
  filter(p1_id == "ieu-a-806", summary_rg_p < 0.05) %>%
  pull(p2_id) %>%
  unique()

selected_traits <- ldsc_results_annotated %>%
  filter(p1_id %in% psychiatric_id, summary_rg_p < 0.05) %>% 
  # filter(p1_id == "ieu-a-806", summary_rg_p < 0.05) %>%
  pull(p2_id) %>%
  unique()


selected_traits <- unique(c(psychiatric_id, selected_traits))

selected_traits <- selected_traits[c(1:20)]

# create_rg_heatmap(
#   selected_traits = selected_traits,
#   ldsc_results_annotated = ldsc_results_annotated,
#   metadata_ieu_EURsampleSize10000 = metadata_ieu_EURsampleSize10000,
#   absolute_correlation = TRUE,
#   text_white_high_threshold = 0.9,
#   text_white_low_threshold = -0.9,
#   rg_text_y_offset_mm = 1.0,
#   star_text_y_offset_mm = -2.5,
#   label_fields_rows = c("id", "trait", "subcategory", "nsnp", "sample_size"),
#   label_fields_cols = c("id", "trait"),
# 
#   row_dend_width_mm = 28,
#   column_dend_height_mm = 24,
# 
#   row_names_max_width_mm_cap = 185,
#   column_names_max_height_mm_cap = 170,
#   
#   output_svg = "/home/mateusz/projects/ifpan-michkor-gr/tmp/bigHeatmap_ldsc_traits_27.03.2026.svg",
#   svg_width = 40,
#   svg_height = 40
# )

create_rg_heatmap(
  selected_traits = selected_traits[c(1:30)],
  ldsc_results_annotated = ldsc_results_annotated,
  metadata_ieu_EURsampleSize10000 = metadata_ieu_EURsampleSize10000,
  absolute_correlation = TRUE,
  text_white_high_threshold = 0.9,
  text_white_low_threshold = -0.9,
  rg_text_y_offset_mm = 1.0,
  star_text_y_offset_mm = -1.5,
  rg_size = 6.5,
  star_size = 8,
  label_fields_rows = c("id", "trait", "subcategory", "nsnp", "sample_size"),
  label_fields_cols = c("id", "trait"),
  
  row_dend_width_mm = 28,
  column_dend_height_mm = 24,
  show_rg_values = F,
  show_significance_stars = F,
  
  row_names_max_width_mm_cap = 185,
  column_names_max_height_mm_cap = 170,
  
  # output_png = "/home/mateusz/projects/ifpan-michkor-gr/tmp/bigHeatmap_ldsc_traits_27.03.2026.png",
  # png_width = c(100, "in"),
  # png_height = c(80, "in"),
  # png_res = 300,
  legend_gap_mm = 20,
  distance_method = "euclidean",
  clustering_method = "average"
)



distance_methods <- c("euclidean", "maximum", "manhattan", "canberra", "binary")
clustering_methods <- c("complete", "average", "single", "ward.D2")

for (dist_method in distance_methods) {
  for (clust_method in clustering_methods) {
    
    output_file <- paste0(
      "/home/mateusz/projects/ifpan-michkor-gr/tmp/bigHeatmap_",
      dist_method, "_", clust_method, ".png"
    )
    
    message("Running: ", dist_method, " + ", clust_method)
    
    create_rg_heatmap(
      selected_traits = selected_traits,
      ldsc_results_annotated = ldsc_results_annotated,
      metadata_ieu_EURsampleSize10000 = metadata_ieu_EURsampleSize10000,
      absolute_correlation = TRUE,
      text_white_high_threshold = 0.9,
      text_white_low_threshold = -0.9,
      rg_text_y_offset_mm = 1.0,
      star_text_y_offset_mm = -1.5,
      rg_size = 6.5,
      star_size = 8,
      label_fields_rows = c("id", "trait", "subcategory", "nsnp", "sample_size"),
      label_fields_cols = c("id", "trait"),
      row_dend_width_mm = 28,
      column_dend_height_mm = 24,
      row_names_max_width_mm_cap = 185,
      column_names_max_height_mm_cap = 170,
      output_png = output_file,
      png_width = c(100, "in"),
      png_height = c(80, "in"),
      png_res = 300,
      legend_gap_mm = 20,
      distance_method = dist_method,
      clustering_method = clust_method
    )
    
  }
}

create_rg_heatmap(
  selected_traits = selected_traits,
  ldsc_results_annotated = ldsc_results_annotated,
  metadata_ieu_EURsampleSize10000 = metadata_ieu_EURsampleSize10000,
  absolute_correlation = TRUE,
  text_white_high_threshold = 0.9,
  text_white_low_threshold = -0.9,
  rg_text_y_offset_mm = 1.0,
  star_text_y_offset_mm = -1.5,
  rg_size = 6.5,
  star_size = 8,
  label_fields_rows = c("id", "trait", "subcategory", "nsnp", "sample_size"),
  label_fields_cols = c("id", "trait"),
  row_dend_width_mm = 28,
  column_dend_height_mm = 24,
  row_names_max_width_mm_cap = 185,
  column_names_max_height_mm_cap = 170,
  output_png = output_file,
  png_width = c(100, "in"),
  png_height = c(80, "in"),
  png_res = 300,
  legend_gap_mm = 20,
  distance_method = dist_method,
  clustering_method = clust_method
)

