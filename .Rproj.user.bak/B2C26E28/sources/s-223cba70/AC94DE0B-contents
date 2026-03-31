################################################################################
# script resposible for prepare visualization enrichr analysis clusters vs chea
################################################################################
chea_list <-
  enrichr_multiple_list_and_databases(data = gr_database_blocked_gene_lists$marpiech_cluster_dex,
                                      databases = c("ChEA_2022"))

chea_list %>% 
  unlist(recursive = F) %>% 
  bind_rows(., .id = "cluster") %>% 
  mutate(cluster = str_replace(cluster, ".ChEA_2022", "")) %>% 
  mutate(cluster = case_when(
    cluster == "cluster_1" ~ "cluster_A",
    cluster == "cluster_2" ~ "cluster_B",
    cluster == "cluster_3" ~ "cluster_C",
    cluster == "cluster_4" ~ "cluster_D",
    cluster == "cluster_5" ~ "cluster_E",
    cluster == "cluster_6" ~ "cluster_F",
    cluster == "cluster_7" ~ "cluster_G",
    cluster == "cluster_8" ~ "cluster_H",
    cluster == "cluster_9" ~ "cluster_I",
    cluster == "cluster_10" ~ "cluster_J",
    cluster == "cluster_11" ~ "cluster_K",
    cluster == "cluster_12" ~ "cluster_L",
    cluster == "cluster_13" ~ "cluster_M",
    cluster == "cluster_14" ~ "cluster_N",
    cluster == "cluster_15" ~ "cluster_O",
    cluster == "cluster_16" ~ "cluster_P",
    cluster == "cluster_17" ~ "DOWN",
    cluster == "cluster_18" ~ "UP",
  )) -> chea_df_raw

chea_df_preprocessing %>% 
  filter(p_value <  0.05) %>% 
  filter(cluster == "cluster_B")


# chea_df_raw %>% 
#   # filter(grepl("STAT3|PPARG|PPARA|RXR|NR3C1", Term, ignore.case = T)) %>% 
#   filter(grepl("RXR|LXR|PPARA|RELB|NR3C1", Term, ignore.case = T)) %>% 
#   mutate(n_genes = str_split_fixed(Overlap, "/", 2)[, 1]) %>% 
#   filter(n_genes > 1) %>% 
#   rename(p_value = "P.value") %>% 
#   rename(fdr = "Adjusted.P.value") %>% 
#   select(c(cluster, Term, p_value, fdr, n_genes ))  %>% 
#   # filter(grepl("STAT3|PPARG|PPARA|RXR|NR3C1", Term, ignore.case = T)) %>% 
#   filter(grepl("RXR|LXR|PPARA|RELB|NR3C1", Term, ignore.case = T)) %>% 
#   # filter(cluster == "cluster_M") %>% 
#   mutate(TF = str_split_fixed(Term, " ", 2)[, 1]) %>% 
#   filter(TF %in% c("RXR", "LXR", "PPARA", "RELB", "NR3C1")) %>% 
#   # filter(p_value < 0.05) %>%
#   mutate(category = ifelse(p_value < 0.05, 1, 0)) %>% 
#   mutate(category = ifelse(fdr < 0.1, 2, category)) %>% 
#   mutate(category = ifelse(fdr < 0.01, 3, category)) %>%
#   group_by(cluster, TF) %>% 
#   filter(p_value == min(p_value)) -> tmp
# 
# tmp %>% 
#   select(-c(Term, p_value, fdr, n_genes)) %>% 
#   spread(key = TF, value = category) %>%
#   replace_na(list(
#     RXR = 0,
#     RELB = 0,
#     PPARA = 0,
#     LXR = 0,
#     NR3C1 = 0
#   )) %>% 
#   melt(., .id = "TF") -> to_plot




create_heatmap <-
  function(data,
           x_axis,
           y_axis,
           palette,
           mapper,
           legend_description,
           x_order = NULL, # New argument for column order
           panel_spacing_size = 6,
           border_size = 1.5,
           tile_size = 0.7,
           border_color = "black",
           x_axis_text_size = 12,
           y_axis_text_size = 12) {
    
    # Sort y-axis elements in reverse alphabetical order and adjust as a factor
    data <- data %>%
      arrange(desc(!!sym(y_axis))) %>%
      mutate(!!y_axis := factor(!!sym(y_axis), levels = unique(!!sym(y_axis))))
    
    # Apply the mapper to create a new column for labeling
    data <- data %>%
      mutate(mapped_value = factor(value, levels = names(mapper), labels = mapper))
    
    # If x_order is provided, adjust the factor levels of x_axis
    if (!is.null(x_order)) {
      data <- data %>%
        mutate(!!x_axis := factor(!!sym(x_axis), levels = x_order))
    }
    
    # Create the heatmap
    plot <-
      ggplot(data, aes(
        x = get(x_axis),
        y = get(y_axis),
        fill = mapped_value
      )) +
      geom_tile(
        color = border_color,
        size = border_size,
        width = tile_size,
        height = tile_size,
        na.rm = TRUE
      ) + # Tiles with custom sizes
      scale_fill_manual(values = palette) + # Use the custom palette passed as argument
      labs(fill = legend_description) +  # Add custom legend title from argument
      theme_minimal() + # Minimal theme
      theme(
        axis.text.x.top = element_text(
          angle = 90,
          hjust = 1,
          size = x_axis_text_size,
          color = "black",
          vjust = 0.5
        ),
        # X-axis labels at the top
        axis.ticks.x.top = element_line(color = "black"),
        # Ensure ticks are at top if needed
        axis.text.x = element_blank(),
        # Hide the default x-axis text at bottom
        axis.text.y = element_text(size = y_axis_text_size),
        # Y-axis text size
        axis.title.x.top = element_blank(),
        # Ensure no title at the top x-axis
        axis.title.x = element_blank(),
        # Hide the bottom x-axis title if it exists
        axis.title.y = element_blank(),
        # Hide y-axis title
        axis.ticks.length = unit(0, "points"),
        # Remove tick marks
        panel.spacing = unit(panel_spacing_size, "lines"),
        # Custom spacing between cells
        plot.background = element_rect(fill = "white", colour = NA),
        # Transparent plot background
        panel.grid.major = element_blank(),
        # Remove major grid lines
        panel.grid.minor = element_blank(),
        # Remove minor grid lines
        panel.border = element_rect(color = "black", fill = NA, size = 1), # Add black border
        legend.position = "bottom",
        # Move legend to bottom
        legend.box = "horizontal",
        # Align legend items horizontally
        legend.title = element_blank() # Adjust legend title size
      ) +
      scale_x_discrete(position = "top")  # Move the x-axis labels to the top
    
    # Add a frame around the entire heatmap
    n_x <- length(unique(data[[x_axis]]))  # Number of x-axis elements
    n_y <- length(unique(data[[y_axis]]))  # Number of y-axis elements
    
    plot
    
  }



################################################################################
# text plot
chea_df_raw %>%
  mutate(n_genes = str_split_fixed(Overlap, "/", 2)[, 1]) %>%
  mutate(n_genes = as.numeric(n_genes)) %>% 
  filter(n_genes >= 2) %>% 
  rename(p_value = "P.value") %>% 
  rename(fdr = "Adjusted.P.value") %>% 
  filter(p_value < 0.05) %>% 
  select(c(cluster, Term, p_value, fdr, n_genes, Genes ))  %>% 
  # mutate(Term = str_replace(Term, "CTB1 Human Placenta Inflammation", "CTB1_(Placenta_Inflammation) Human")) %>% 
  # mutate(Term = str_replace(Term, "C42B Human ProstateCancer", "C42B_(Prostate_Cancer) Human")) %>% 
  # mutate(Term = str_replace(Term, "H9 Human BoneMarrow Lymphoma", "H9_(Bone_Marrow_Lymphoma) Human")) %>% 
  # mutate(Term = str_replace(Term, "WistarRat Hippocampus Stress", "Hippocampus_(Stress) WistarRat")) %>% 
  # mutate(Term = str_replace(Term, "WistarRat Hippocampus", "Hippocampus WistarRat")) %>% 
  # mutate(Term = str_replace(Term, "A549 Human Lung Carcinoma", "A549_(Lung_Carcinoma) Human")) %>% 
  # mutate(Term = str_replace(Term, "MC3T3E1 Mouse Bone", "MC3T3E1_(Bone) Mouse")) %>% 
  mutate(TF = str_split_fixed(Term, " ", 5)[, 1]) %>% 
  mutate(PMID = str_split_fixed(Term, " ", 5)[, 2]) %>% 
  mutate(method = str_split_fixed(Term, " ", 5)[, 3]) %>% 
  mutate(tissue = str_split_fixed(Term, " ", 5)[, 4]) %>% 
  mutate(species = str_split_fixed(Term, " ", 5)[, 5]) %>% 
  group_by(cluster, TF) %>% 
  slice_min(p_value) %>%
  ungroup %>% 
  group_by(cluster) %>% 
  arrange(p_value) %>% 
  slice_head(n = 3) %>% 
  mutate(top = row_number()) %>% 
  select(c(cluster, Term, top, p_value, fdr, TF, PMID, tissue, species)) %>% 
  mutate(label = paste(top, TF, sep = ". ")) %>% 
  as.data.frame()  %>% 
  mutate(id = rep(c(1,2.5,4),18)) %>% 
  mutate(cluster_letters = str_replace_all(cluster, "cluster_", "")) -> to_text_plot

to_text_plot %>%
  mutate(text_properties = ifelse(TF %in% c("NR3C1","RXR", "LXR", "PPARA", "RELB"), "bold", "plain")) %>%
  mutate(significant = ifelse(p_value < 0.05, "*", "n.s")) %>% 
  mutate(significant = ifelse(fdr < 0.1, "**", significant)) %>% 
  mutate(significant = ifelse(fdr < 0.01, "***", significant)) %>% 
  mutate(label2 = paste(label, significant, sep = "  ")) -> to_text_plot

ggplot() +
  # geom_text(data = to_text_plot, aes(x = 1.0, y = id, label = TF), hjust = 0) +
  geom_text(data = filter(to_text_plot, text_properties == "plain"), aes(x = 1.0, y = id, label = label), hjust = 0, family = "sans") +
  geom_text(data = filter(to_text_plot, text_properties == "bold"), aes(x = 1.0, y = id, label = label), hjust = 0, fontface = "bold", family = "sans") +
  facet_grid(cluster ~ .) +  # Facet grid for different clusters
  scale_x_continuous(limits = c(1, 1.02)) +  # Set y-axis limits
  scale_y_reverse(limits = c(5, 0)) + 
  theme_classic() +
  theme(
    axis.text.x = element_blank(),  # Hide x-axis text
    axis.ticks.x = element_blank(),  # Hide x-axis ticks
    axis.title.x = element_blank(),  # Hide x-axis title
    axis.title.y = element_blank(),  # Hide y-axis title
    axis.text.y = element_blank(),   # Hide y-axis text
    axis.ticks.y = element_blank(),  # Hide y-axis ticks
    strip.text = element_blank()     # Hide facet titles
  ) -> text_plot

svg("results/google-drive/enrichr/gr-dependent-transcriptional-pattern/chea_text_plot.svg", width = 5, height = 13)
text_plot
dev.off()



# chea plot version 2

chea_df_raw %>% 
  # filter(grepl("STAT3|PPARG|PPARA|RXR|NR3C1", Term, ignore.case = T)) %>% 
  filter(grepl("RXR|LXR|PPARA|RELB|NR3C1|SUZ12|IRF8", Term, ignore.case = T)) %>% 
  mutate(n_genes = str_split_fixed(Overlap, "/", 2)[, 1]) %>% 
  filter(n_genes > 1) %>% 
  rename(p_value = "P.value") %>% 
  rename(fdr = "Adjusted.P.value") %>% 
  select(c(cluster, Term, p_value, fdr, n_genes ))  %>% 
  # filter(grepl("STAT3|PPARG|PPARA|RXR|NR3C1", Term, ignore.case = T)) %>% 
  filter(grepl("RXR|LXR|PPARA|RELB|NR3C1|SUZ12|IRF8", Term, ignore.case = T)) %>% 
  # filter(cluster == "cluster_M") %>% 
  mutate(TF = str_split_fixed(Term, " ", 2)[, 1]) %>% 
  filter(TF %in% c("RXR", "LXR", "PPARA", "RELB", "NR3C1", "IRF8", "SUZ12")) %>% 
  # filter(p_value < 0.05) %>%
  mutate(category = ifelse(p_value < 0.05, 1, 0)) %>% 
  mutate(category = ifelse(fdr < 0.1, 2, category)) %>% 
  mutate(category = ifelse(fdr < 0.01, 3, category)) %>%
  group_by(cluster, TF) %>% 
  filter(p_value == min(p_value)) %>% 
  select(-c(Term, p_value, fdr, n_genes)) %>% 
  spread(key = TF, value = category) %>%
  replace_na(list(
    RXR = 0,
    RELB = 0,
    PPARA = 0,
    LXR = 0,
    NR3C1 = 0,
    IRF8 = 0,
    SUZ12 = 0
  )) %>% 
  melt(., .id = "TF") -> to_plot2



create_heatmap(
  data = to_plot2,
  x_axis = "variable",
  y_axis = "cluster",
  palette = palette,
  border_color = "white",
  mapper = mapper_legend,
  legend_description = "",
  border_size = 1.5,
  tile_size = 1,
  x_order = c("NR3C1", "RXR", "LXR", "PPARA", "RELB", "SUZ12", "IRF8"),
  x_axis_text_size = 14,
  y_axis_text_size = 14
)


# palette <- c("0" = "#e2e2e2", "1" = "#658354", "2" = "#e6352b", "3" = "#e89796ff")
mapper_legend <-  c("0" = "n.s.", "1" = "p < 0.05", "2" = "FDR < 0.1", "3" = "FDR < 0.01")
palette <- c("n.s." = "#ffffffff", "p < 0.05" = "#fce3e2ff", "FDR < 0.1" = "#f4bcbaff", "FDR < 0.01" = "#e68b8aff")

svg("results/google-drive/enrichr/gr-dependent-transcriptional-pattern/chea_TF_heatmap.svg", width = 6, height = 18)

to_plot %>% 
  select(cluster) %>% 
  unique() %>% 
  mutate(label = c("cluster A (9)", "cluster B (7)", "cluster C (19)", "cluster D (40)",
                   "cluster E (7)", "cluster F (11)", "cluster G (14)", "cluster H (15)",
                   "cluster I (27)", "cluster J (18)", "cluster K (25)", "cluster L (75)",
                   "cluster M (13)", "cluster N (6)", "cluster O (54)", "cluster P (21)",
                   "DOWN (93)", "UP (269)")) %>% 
  deframe -> mapper_cluster

to_plot %>% 
  left_join(., mapper_label_cluster, by = "cluster") -> to_plot

create_heatmap(
  data = to_plot,
  x_axis = "variable",
  y_axis = "cluster",
  palette = palette,
  border_color = "white",
  legend_mapper = mapper_legend,
  legend_description = "",
  border_size = 1.5,
  tile_size = 1,
  x_order = c("NR3C1", "RXR", "LXR", "PPARA", "RELB"),
  x_axis_text_size = 14,
  y_axis_text_size = 14
)

dev.off()


################################################################################
# rotate heatmap for TF
create_heatmap <- function(data,
                           x_axis,
                           y_axis,
                           palette,
                           legend_mapper,
                           legend_description,
                           x_order = NULL,
                           y_order = NULL,
                           x_axis_mapper = NULL,
                           x_axis_rotate = 45, # Obrót tekstu dla lepszej czytelności
                           panel_spacing_size = 6,
                           border_size = 1.5,
                           tile_size = 0.7,
                           border_color = "black",
                           x_axis_text_size = 12, # Dopasowany rozmiar tekstu
                           y_axis_text_size = 12) {
  
  # Sortowanie osi Y
  data <- data %>%
    arrange(desc(!!sym(y_axis))) %>%
    mutate(!!y_axis := factor(!!sym(y_axis), levels = unique(!!sym(y_axis))))
  
  # Mapa legendy
  data <- data %>%
    mutate(mapped_value = factor(value, levels = names(legend_mapper), labels = legend_mapper))
  
  # Ustawienie kolejności osi X
  if (!is.null(x_order)) {
    data <- data %>%
      mutate(!!x_axis := factor(!!sym(x_axis), levels = x_order))
  }
  
  # Ustawienie kolejności osi Y
  if (!is.null(y_order)) {
    data <- data %>%
      mutate(!!y_axis := factor(!!sym(y_axis), levels = y_order))
  }
  
  # Zmapowanie tekstu osi X
  if (!is.null(x_axis_mapper)) {
    data <- data %>%
      mutate(!!x_axis := factor(!!sym(x_axis), levels = names(x_axis_mapper), labels = x_axis_mapper))
  }
  
  # Tworzenie heatmapy
  plot <- ggplot(data, aes(x = get(x_axis), y = get(y_axis), fill = mapped_value)) +
    geom_tile(color = border_color, size = border_size, width = tile_size, height = tile_size, na.rm = TRUE) +
    scale_fill_manual(values = palette) +
    labs(fill = legend_description) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = x_axis_rotate, size = x_axis_text_size, vjust = 0, hjust = 0), # Ustawienie na górze z kątem
      axis.text.y = element_text(size = y_axis_text_size),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.length = unit(0, "points"),
      panel.spacing = unit(panel_spacing_size, "lines"),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_blank()
    ) +
    scale_x_discrete(position = "top") # Tekst osi X na górze
  
  plot
}


################################################################################
svg("results/google-drive/enrichr/gr-dependent-transcriptional-pattern/chea_TF_heatmap-horizontal.svg", width = 18, height = 7.4)

create_heatmap(
  data = to_plot %>% mutate(label = "test"),
  x_axis = "cluster",
  y_axis = "variable",
  palette = palette,
  border_color = "white",
  legend_mapper = mapper_legend,
  x_axis_mapper = mapper_cluster,
  x_axis_rotate = 45,
  legend_description = "",
  border_size = 1.5,
  tile_size = 1,
  y_order = rev(c("NR3C1", "RXR", "LXR", "PPARA", "RELB")),
  x_axis_text_size = 28,
  y_axis_text_size = 28
) 

dev.off()

to_text_plot %>% head

label_size = 8
ggplot() +
  geom_text(data = filter(to_text_plot, text_properties == "plain"), 
            aes(x = id, y = 1.0, label = TF), 
            hjust = 0, 
            family = "sans", size = label_size) +
  geom_text(data = filter(to_text_plot, text_properties == "bold"), 
            aes(x = id, y = 1.0, label = TF), 
            hjust = 0, 
            fontface = "bold", 
            family = "sans",
            size = label_size) +
  facet_grid(. ~ cluster) +  # Adjust facet layout to horizontal
  scale_x_reverse(limits = c(5, 0)) +  # Reverse x-axis for horizontal layout
  scale_y_continuous(limits = c(1, 1.02)) +  # Set y-axis limits for spacing
  theme_classic() +
  theme(
    axis.text.x = element_blank(),  # Hide x-axis text
    axis.ticks.x = element_blank(),  # Hide x-axis ticks
    axis.title.x = element_blank(),  # Hide x-axis title
    axis.title.y = element_blank(),  # Hide y-axis title
    axis.text.y = element_blank(),   # Hide y-axis text
    axis.ticks.y = element_blank(),  # Hide y-axis ticks
    strip.text = element_blank()     # Hide facet titles
  ) +
  coord_flip() -> text_plot

svg("results/google-drive/enrichr/gr-dependent-transcriptional-pattern/chea_text_plot-horizontal.svg", width = 22, height = 2)
text_plot
dev.off()



###################################################################
library(ComplexHeatmap)

# Przykładowa heatmapa
mat <- matrix(rnorm(100), 10, 10)
ht <- Heatmap(mat)

# Sprawdź bieżące ustawienia dotyczące czcionki dla nazw wierszy i kolumn
ht@row_names_param$gp
ht@column_names_param$gp

# Sprawdź bieżący motyw
current_theme <- theme_get()

# Wyświetl elementy tekstowe z motywu
current_theme$text

