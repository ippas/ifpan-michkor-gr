metadata_ieu_EURsampleSize10000 %>% 
  filter(subcategory == "Psychiatric / neurological") %>% 
  filter(trait %in% c(
    "Major depressive disorder",
    "Major Depressive Disorder",
    "Depressive symptoms",
    "MZ twin differences on depression symptoms score",
    "MZ twin differences on depression symptoms score in children",
    "MZ twin differences on depression symptoms score in adults",
    "Bipolar disorder",
    "bipolar disorder",
    "Bipolar disorder bip2021",
    "Subjective well being",
    "MZ twin differences on subjective wellbeing score",
    "Neuroticism"
  )) %>% .$id -> psychiatric_id




ldsc_results_annotated %>%
  filter(p2_subcategory == "Psychiatric / neurological") %>%
  group_by(p2_trait) %>%
  slice_max(order_by = abs(summary_rg), n = 1, with_ties = FALSE) %>%
  ungroup() %>% 
  pull(p2_id) %>%
  unique() -> psychiatric_id

selected_traits <- ldsc_results_annotated %>%
  filter(p1_id %in% psychiatric_id, summary_rg_p < 0.01) %>% 
  filter(p2_subcategory != "Psychiatric / neurological") %>%
  group_by(p2_trait) %>%
  slice_max(order_by = abs(summary_rg), n = 1, with_ties = FALSE) %>%
  ungroup() %>% 
  # filter(subcategory == "Psychiatric / neurological")
  # filter(p1_id == "ieu-a-806", summary_rg_p < 0.05) %>%
  pull(p2_id) %>%
  unique()

# ldsc_results_annotated %>%
#   filter(p1_id %in% psychiatric_id, summary_rg_p < 0.001) %>% 
#   filter(p2_subcategory != "Psychiatric / neurological") %>% .$p2_trait %>% unique()

# selected_traits <- unique(c(psychiatric_id, selected_traits))

create_rg_heatmap(
  selected_traits = selected_traits,
  ldsc_results_annotated = ldsc_results_annotated,
  metadata_ieu_EURsampleSize10000 = metadata_ieu_EURsampleSize10000,
  absolute_correlation = TRUE,
  text_white_high_threshold = 0.9,
  text_white_low_threshold = -0.9,
  rg_text_y_offset_mm = 1.0,
  star_text_y_offset_mm = -2.5,
  label_fields_rows = c("id", "trait", "subcategory", "nsnp", "sample_size"),
  label_fields_cols = c("id", "trait"),
  
  row_dend_width_mm = 28,
  column_dend_height_mm = 24,
  
  row_names_max_width_mm_cap = 185,
  column_names_max_height_mm_cap = 170,
  show_rg_values = T,
  show_significance_stars = T,
  # color_row_dend_branches = TRUE,
  # color_col_dend_branches = TRUE,
  # row_dend_k = 7,
  # col_dend_k = 7,
  # row_dend_colors = c("#D73027", "#4575B4", "#1A9850", "#984EA3"),
  # col_dend_colors = c("#D73027", "#4575B4", "#1A9850", "#984EA3"),
  color_row_dend_branches = TRUE,
  color_col_dend_branches = TRUE,
  
  row_dend_k = 5,
  col_dend_k = 5,
  row_gap_color = "black",
  col_gap_color = "white",
  use_row_cluster_gaps = TRUE,
  use_col_cluster_gaps = T,
  
  row_gap_width_mm = 1,
  col_gap_width_mm = 1,
  distance_method = "euclidean",
  clustering_method = "complete"
  # output_png = "/home/mateusz/projects/ifpan-michkor-gr/results_v2/ldsc_results/heatmaps/phenotypesAssociatedDepressionTraits/HeatmapLDSCAssociatedDepressionTraits_29.03.2026.png",
  # png_width = c(32, "in"),
  # png_height = c(25, "in"),
  # png_res = 300
)


create_rg_heatmap(
  selected_traits = selected_traits,
  ldsc_results_annotated = ldsc_results_annotated,
  metadata_ieu_EURsampleSize10000 = metadata_ieu_EURsampleSize10000,
  absolute_correlation = TRUE,
  text_white_high_threshold = 0.9,
  text_white_low_threshold = -0.9,
  rg_text_y_offset_mm = 1.0,
  star_text_y_offset_mm = -2.5,
  label_fields_rows = c("id", "trait", "subcategory", "nsnp", "sample_size"),
  label_fields_cols = c("id", "trait"),

  row_dend_width_mm = 28,
  column_dend_height_mm = 24,

  row_names_max_width_mm_cap = 185,
  column_names_max_height_mm_cap = 170,
  show_rg_values = T,
  show_significance_stars = T,
  # color_row_dend_branches = TRUE,
  # color_col_dend_branches = TRUE,
  # row_dend_k = 7,
  # col_dend_k = 7,
  # row_dend_colors = c("#D73027", "#4575B4", "#1A9850", "#984EA3"),
  # col_dend_colors = c("#D73027", "#4575B4", "#1A9850", "#984EA3"),
  color_row_dend_branches = TRUE,
  color_col_dend_branches = TRUE,
  
  row_dend_k = 7,
  col_dend_k = 7,
  row_gap_color = "black",
  col_gap_color = "white",
  use_row_cluster_gaps = TRUE,
  use_col_cluster_gaps = T,
  
  row_gap_width_mm = 1,
  col_gap_width_mm = 1,
  distance_method = "canberra",
  clustering_method = "ward.D2",

  output_png = "/home/mateusz/projects/ifpan-michkor-gr/results_v2/ldsc_results/heatmaps/phenotypesAssociatedDepressionTraits/HeatmapLDSCAssociatedDepressionTraits_29.03.2026.png",
  png_width = c(32, "in"),
  png_height = c(25, "in"),
  png_res = 300
)


distance_methods <- c(
  "euclidean",
  "maximum",
  "manhattan",
  "canberra",
  "binary",
  "minkowski"
)

clustering_methods <- c(
  "complete",
  "single",
  "average",
  "mcquitty",
  "median",
  "centroid",
  "ward.D",
  "ward.D2"
)

cluster_k_values <- 4:10

output_dir <- "/home/mateusz/projects/ifpan-michkor-gr/results_v2/ldsc_results/heatmaps/phenotypesAssociatedDepressionTraits"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

for (distance_method_i in distance_methods) {
  for (clustering_method_i in clustering_methods) {
    for (k_i in cluster_k_values) {
      
      output_file <- file.path(
        output_dir,
        paste0(
          "HeatmapLDSCAssociatedDepressionTraits",
          "_p0.01",
          "_k", k_i,
          "_dist-", distance_method_i,
          "_clust-", clustering_method_i,
          "_29.03.2026.png"
        )
      )
      
      message("========================================")
      message("Running heatmap:")
      message("distance_method = ", distance_method_i)
      message("clustering_method = ", clustering_method_i)
      message("k = ", k_i)
      message("output = ", output_file)
      
      tryCatch(
        {
          create_rg_heatmap(
            selected_traits = selected_traits,
            ldsc_results_annotated = ldsc_results_annotated,
            metadata_ieu_EURsampleSize10000 = metadata_ieu_EURsampleSize10000,
            
            absolute_correlation = TRUE,
            text_white_high_threshold = 0.9,
            text_white_low_threshold = -0.9,
            rg_text_y_offset_mm = 1.0,
            star_text_y_offset_mm = -2.5,
            
            label_fields_rows = c("id", "trait", "subcategory", "nsnp", "sample_size"),
            label_fields_cols = c("id", "trait"),
            
            row_dend_width_mm = 28,
            column_dend_height_mm = 24,
            
            row_names_max_width_mm_cap = 185,
            column_names_max_height_mm_cap = 170,
            
            show_rg_values = TRUE,
            show_significance_stars = TRUE,
            
            color_row_dend_branches = TRUE,
            color_col_dend_branches = TRUE,
            
            row_dend_k = k_i,
            col_dend_k = k_i,
            
            row_gap_color = "black",
            col_gap_color = "white",
            
            use_row_cluster_gaps = TRUE,
            use_col_cluster_gaps = TRUE,
            
            row_gap_width_mm = 1,
            col_gap_width_mm = 1,
            
            distance_method = distance_method_i,
            clustering_method = clustering_method_i,
            
            output_png = output_file,
            png_width = c(32, "in"),
            png_height = c(25, "in"),
            png_res = 300
          )
          
          message("Finished successfully: ", output_file)
        },
        error = function(e) {
          message("ERROR for:")
          message("  distance_method = ", distance_method_i)
          message("  clustering_method = ", clustering_method_i)
          message("  k = ", k_i)
          message("  message = ", conditionMessage(e))
        }
      )
    }
  }
}
