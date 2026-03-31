plot_signature_density <- function(data,
                                   background_list,
                                   signature_names,
                                   background_colors = c("#7F7F7F", "#1B9E77"),
                                   other_colors = NULL,
                                   background_line_size = 1.3,
                                   other_line_size = 0.7,
                                   log_y_axis = TRUE) {
  stopifnot(length(background_list) == length(background_colors))
  
  # Filtrowanie danych
  df <- data %>%
    mutate(signature_derivation = signature_name) %>%
    filter(signature_name %in% signature_names) %>%
    filter(!signature_derivation %in% c("cluster_M")) %>%
    mutate(log_pvalue = -log10(pvalue))
  
  # Gęstości
  density_df <- df %>%
    group_by(signature_derivation) %>%
    summarise(density_data = list({
      d <- density(log_pvalue, adjust = 1.0)
      data.frame(x = d$x, y = d$y)
    }), .groups = "drop") %>%
    unnest(density_data) %>%
    mutate(
      line_size = ifelse(signature_derivation %in% background_list,
                         background_line_size, other_line_size),
      y_value = if (log_y_axis) log10(y + 1e-10) else y
    )
  
  # Kolory
  custom_colors <- setNames(background_colors, background_list)
  all_signatures <- unique(density_df$signature_derivation)
  missing_signatures <- setdiff(all_signatures, names(custom_colors))
  
  if (is.null(other_colors)) {
    if (length(missing_signatures) > 0) {
      auto_colors <- hue_pal()(length(missing_signatures))
      names(auto_colors) <- missing_signatures
    } else {
      auto_colors <- character(0)
    }
  } else {
    auto_colors <- other_colors
    missing_auto <- setdiff(missing_signatures, names(auto_colors))
    if (length(missing_auto) > 0) {
      warning("Brakuje kolorów dla: ", paste(missing_auto, collapse = ", "))
    }
  }
  
  final_colors <- c(custom_colors, auto_colors)
  
  # Etykieta osi Y
  y_label <- if (log_y_axis) expression(log[10]("Gęstość")) else "Gęstość"
  
  # Wykres
  ggplot(density_df, aes(x = x, y = y_value,
                         color = signature_derivation,
                         group = signature_derivation)) +
    geom_line(aes(size = line_size)) +
    scale_color_manual(values = final_colors) +
    scale_size_identity() +
    labs(
      x = expression(-log[10](p)),
      y = y_label,
      color = "Signature Derivation"
    ) +
    theme_minimal()
}

# Wstępne filtrowanie danych
data_filtered <- rbind(
  pgc_mdd2025_no23andMe_eur_processing,
  intesection_results_all_gr_signatures,
  enrichr_depression_genes_preprocessing
) %>% 
  filter(range_plus_minus %in% c("100000", "all"),
         pvalue < 1e-8)

# Stałe
background_list <- c("GWAS", "Mental_Depression_546/575_DisGeNET")
signature_names <- c(background_list, "cluster_A", "cluster_B", "cluster_C", "cluster_D", "cluster_E", "cluster_F")

# Rysowanie
plot_signature_density(data_filtered, background_list, signature_names, log_y_axis = FALSE)


signature_names <- c(background_list, "cluster_G", "cluster_H", "cluster_I", "cluster_J", "cluster_K", "cluster_L", "cluster_M", "cluster_N", "cluster_O", "cluster_P")

# Rysowanie
plot_signature_density(data_filtered, background_list, signature_names)



################################################################################
signature_names <- c(background_list, "universal_up", "universal_down")

other_colors <- c(
  "universal_up" = "firebrick",
  "universal_down" = "blue4"
)

plot_signature_density(
  data = data_filtered,
  background_list = background_list,
  signature_names = signature_names,
  background_colors = c("#7F7F7F", "#1B9E77"),  # np. GWAS, Mental_Depression...
  other_colors = other_colors,
  log_y_axis = FALSE
)

################################################################################
signature_names <- c(background_list,
                     "microglia_up", "embryos_up", "lymphoid_up", "muscle_up",
                     "macrophage-like_up", "lung_up", "adrenal-gland_up",
                     "epithelial_up", "ASM_up", "cartilage_up")

# other_colors <- c(
#   "universal_up" = "firebrick",
#   "universal_down" = "blue4"
# )

plot_signature_density(
  data = data_filtered,
  background_list = background_list,
  signature_names = signature_names,
  background_colors = c("#7F7F7F", "#1B9E77"),  # np. GWAS, Mental_Depression...
  # other_colors = other_colors,
  log_y_axis = FALSE
)


signature_names <- c(background_list,
                     "blood_down", "lymphoid_down", "ASM_down", "placenta_down", "lung_down",
                     "bone_down", "epithelial_down", "macrophage-like_down", "spleen_down",
                     "muscle_down", "adrenal-gland_down", "cartilage_down",
                     "kidney_down")

# other_colors <- c(
#   "universal_up" = "firebrick",
#   "universal_down" = "blue4"
# )

plot_signature_density(
  data = data_filtered,
  background_list = background_list,
  signature_names = signature_names,
  background_colors = c("#7F7F7F", "#1B9E77"),  # np. GWAS, Mental_Depression...
  # other_colors = other_colors,
  log_y_axis = FALSE
)
