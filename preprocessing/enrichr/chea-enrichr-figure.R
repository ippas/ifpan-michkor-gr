# 
# create_heatmap <- function(data, x_axis, y_axis, panel_spacing_size = 6, border_size = 1.5, tile_size = 0.7,
#                            x_axis_text_size = 12, y_axis_text_size = 12) {
#   # Adjusting the factor levels to add space between specific columns
#   # data$variable <- factor(data$variable, levels = c("Phenotype", "blank", "Expression", "Localization", "Pathway"))
#   
#   # Create the heatmap
#   plot <- ggplot(data, aes(x = get(x_axis), y = get(y_axis), fill = factor(value))) +
#     geom_tile(color = "black", size = border_size, width = tile_size, height = tile_size, na.rm = TRUE) + # Tiles with custom sizes
#     scale_fill_manual(values = c("0" = "grey", "1" = "orange", "2" = "red")) + # Custom colors
#     theme_minimal() + # Minimal theme
#     theme(
#       axis.text.x.top = element_text(angle = 90, hjust = 1, size = x_axis_text_size, color = "black", vjust = 0.5), # X-axis labels at the top
#       axis.ticks.x.top = element_line(color = "black"), # Ensure ticks are at top if needed
#       axis.text.x = element_blank(), # Hide the default x-axis text at bottom
#       axis.text.y = element_text(size = y_axis_text_size), # Y-axis text size
#       axis.title.x.top = element_blank(), # Ensure no title at the top x-axis
#       axis.title.x = element_blank(), # Hide the bottom x-axis title if it exists
#       axis.title.y = element_blank(), # Hide y-axis title
#       axis.ticks.length = unit(0, "points"), # Remove tick marks
#       panel.spacing = unit(panel_spacing_size, "lines"), # Custom spacing between cells
#       plot.background = element_rect(fill = "white", colour = NA), # Transparent plot background
#       panel.grid.major = element_blank(), # Remove major grid lines
#       panel.grid.minor = element_blank(), # Remove minor grid lines
#       strip.background = element_rect(fill = "white", colour = "black"), # Adjust the background of facet labels
#       panel.border = element_blank(), # Remove default panel border
#       legend.position = "bottom", # Move legend to bottom
#       legend.box = "horizontal", # Align legend items horizontally
#       legend.title = element_blank() # Remove legend title if desired
#     ) +
#     scale_x_discrete(position = "top")  # Move the x-axis labels to the top
#   
#   # Add a vertical line on the left side only
#   plot + geom_vline(xintercept = 0.5, color = "black", size = 1)
# }


enrichr_all_database_gr_dependent_transcriptional_pattern_df %>% 
  filter(database_name %in% c("ChEA_2022")) %>% as.data.frame() %>% 
  filter(n_genes >= 2, P.value < 0.05) %>% 
  filter(database_name == "GO_Biological_Process_2023") 


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

chea_df_raw



chea_df_raw %>%
  filter(grepl("STAT3|PPARG|PPARA|RXR|NR3C1", Term, ignore.case = T)) %>%
  mutate(n_genes = str_split_fixed(Overlap, "/", 2)[, 1]) %>%
  filter(n_genes > 1) %>%
  rename(p_value = "P.value") %>%
  rename(fdr = "Adjusted.P.value") %>%
  select(c(cluster, Term, p_value, fdr, n_genes )) %>%
  mutate(
    contains_STAT3 = if_else(grepl("STAT3", Term, ignore.case = TRUE), 1, 0),
    contains_PPARG = if_else(grepl("PPARG", Term, ignore.case = TRUE), 1, 0),
    contains_PPARA = if_else(grepl("PPARA", Term, ignore.case = TRUE), 1, 0),
    contains_RXR = if_else(grepl("RXR", Term, ignore.case = TRUE), 1, 0),
    contains_NR3C1 = if_else(grepl("NR3C1", Term, ignore.case = TRUE), 1, 0)
  ) -> chea_df_preprocessing


chea_df_raw %>%
  filter(grepl("LXR|RELB|PPARA|RXR|NR3C1", Term, ignore.case = T)) %>%
  mutate(n_genes = str_split_fixed(Overlap, "/", 2)[, 1]) %>%
  filter(n_genes > 1) %>%
  rename(p_value = "P.value") %>%
  rename(fdr = "Adjusted.P.value") %>%
  select(c(cluster, Term, p_value, fdr, n_genes )) %>%
  mutate(
    contains_STAT3 = if_else(grepl("RXR", Term, ignore.case = TRUE), 1, 0),
    contains_PPARG = if_else(grepl("LXR", Term, ignore.case = TRUE), 1, 0),
    contains_PPARA = if_else(grepl("PPARA", Term, ignore.case = TRUE), 1, 0),
    contains_RXR = if_else(grepl("RELB", Term, ignore.case = TRUE), 1, 0),
    contains_NR3C1 = if_else(grepl("NR3C1", Term, ignore.case = TRUE), 1, 0)
  ) -> chea_df_preprocessing


chea_df_preprocessing %>% 
  filter(p_value <  0.05) %>% 
  gather(., key = "interest_TF", value = "contain_TF", -c(cluster, Term, p_value, fdr, n_genes)) %>% 
  mutate(contain_TF = ifelse(fdr < 0.05, 2, contain_TF)) %>% 
  group_by(cluster, interest_TF) %>%
  nest() %>%
  mutate(
    data = map(data, ~ .x %>% arrange(p_value) %>% slice(1))
  ) %>%
  unnest(data) %>% 
  # # filter(!(cluster %in% c("UP", "DOWN"))) %>%
  filter(contain_TF != 0) %>%
  select(-Term) %>%
  ungroup %>%
  group_by(cluster, interest_TF) %>%
  slice_max(contain_TF, with_ties = TRUE) %>% 
  spread(key = interest_TF, value = contain_TF) %>%
  replace_na(list(
    contains_STAT3 = 0,
    contains_PPARG = 0,
    contains_PPARA = 0,
    contains_RXR = 0,
    contains_NR3C1 = 0
  )) %>% 
  select(-c(p_value, fdr)) %>% 
  melt(., .id = "interest_TF") %>% 
  mutate(variable = str_replace(variable, "contains_", "")) -> chea_to_plot



create_heatmap(data=chea_to_plot, x_axis = "variable", y_axis = "cluster",
               border_size = 1, tile_size = 0.7, x_axis_text_size = 14, y_axis_text_size = 14)


create_heatmap <- function(data, x_axis, y_axis, panel_spacing_size = 6, border_size = 1.5, tile_size = 0.7,
                           x_axis_text_size = 12, y_axis_text_size = 12) {
  # Sort y-axis elements in reverse alphabetical order and adjust as a factor
  data <- data %>%
    arrange(desc(!!sym(y_axis))) %>%
    mutate(!!y_axis := factor(!!sym(y_axis), levels = unique(!!sym(y_axis))))
  
  # Create the heatmap
  plot <- ggplot(data, aes(x = get(x_axis), y = get(y_axis), fill = factor(value))) +
    geom_tile(color = "black", size = border_size, width = tile_size, height = tile_size, na.rm = TRUE) + # Tiles with custom sizes
    scale_fill_manual(values = c("0" = "#e2e2e2", "1" = "#ff6242", "2" = "#fb3b1e", "3" = "#c61a09")) + # Custom colors
    # scale_fill_manual(values = c("0" = "grey", "1" = "orange", "2" = "red")) + # Custom colors
    theme_minimal() + # Minimal theme
    theme(
      axis.text.x.top = element_text(angle = 90, hjust = 1, size = x_axis_text_size, color = "black", vjust = 0.5), # X-axis labels at the top
      axis.ticks.x.top = element_line(color = "black"), # Ensure ticks are at top if needed
      axis.text.x = element_blank(), # Hide the default x-axis text at bottom
      axis.text.y = element_text(size = y_axis_text_size), # Y-axis text size
      axis.title.x.top = element_blank(), # Ensure no title at the top x-axis
      axis.title.x = element_blank(), # Hide the bottom x-axis title if it exists
      axis.title.y = element_blank(), # Hide y-axis title
      axis.ticks.length = unit(0, "points"), # Remove tick marks
      panel.spacing = unit(panel_spacing_size, "lines"), # Custom spacing between cells
      plot.background = element_rect(fill = "white", colour = NA), # Transparent plot background
      panel.grid.major = element_blank(), # Remove major grid lines
      panel.grid.minor = element_blank(), # Remove minor grid lines
      strip.background = element_rect(fill = "white", colour = "black"), # Adjust the background of facet labels
      panel.border = element_blank(), # Remove default panel border
      legend.position = "bottom", # Move legend to bottom
      legend.box = "horizontal", # Align legend items horizontally
      legend.title = element_blank() # Remove legend title if desired
    ) +
    scale_x_discrete(position = "top")  # Move the x-axis labels to the top
  
  # Add a vertical line on the left side only
  plot + geom_vline(xintercept = 0.5, color = "black", size = 1)
}



chea_df_preprocessing %>% 
  filter(p_value <  0.05) %>% 
  filter(cluster == "cluster_B")


chea_df_raw %>% 
  filter(grepl("STAT3|PPARG|PPARA|RXR|NR3C1", Term, ignore.case = T)) %>% 
  mutate(n_genes = str_split_fixed(Overlap, "/", 2)[, 1]) %>% 
  filter(n_genes > 1) %>% 
  rename(p_value = "P.value") %>% 
  rename(fdr = "Adjusted.P.value") %>% 
  select(c(cluster, Term, p_value, fdr, n_genes ))  %>% 
  filter(grepl("STAT3|PPARG|PPARA|RXR|NR3C1", Term, ignore.case = T)) %>% 
  # filter(cluster == "cluster_M") %>% 
  mutate(TF = str_split_fixed(Term, " ", 2)[, 1]) %>% 
  filter(TF %in% c("STAT3", "PPARG", "PPARA", "RXR", "NR3C1")) %>% 
  # filter(p_value < 0.05) %>%
  mutate(category = ifelse(p_value < 0.05, 1, 0)) %>% 
  mutate(category = ifelse(fdr < 0.1, 2, category)) %>% 
  mutate(category = ifelse(fdr < 0.01, 3, category)) %>%
  group_by(cluster, TF) %>% 
  filter(p_value == min(p_value)) -> tmp

tmp %>% 
  select(-c(Term, p_value, fdr, n_genes)) %>% 
  spread(key = TF, value = category) %>%
  replace_na(list(
    STAT3 = 0,
    PPARG = 0,
    PPARA = 0,
    RXR = 0,
    NR3C1 = 0
  )) %>% 
  melt(., .id = "TF") -> to_plot

to_plot %>% filter(cluster == "cluster_")

svg("results/google-drive/enrichr/gr-dependent-transcriptional-pattern/chea_TF_heatmap2.svg", width = 5, height = 15)
create_heatmap(data=to_plot, x_axis = "variable", y_axis = "cluster",
               border_size = 1, tile_size = 0.6, x_axis_text_size = 14, y_axis_text_size = 14)

dev.off()


chea_df_raw %>% 
  mutate(n_genes = str_split_fixed(Overlap, "/", 2)[, 1]) %>% 
  filter(n_genes >= 2) %>% 
  rename(p_value = "P.value") %>% 
  rename(fdr = "Adjusted.P.value") %>% 
  filter(p_value < 0.05) %>% 
  select(c(cluster, Term, p_value, fdr, n_genes, Genes ))  %>% 
  mutate(Term = str_replace(Term, "CTB1 Human Placenta Inflammation", "CTB1_(Placenta_Inflammation) Human")) %>% 
  mutate(Term = str_replace(Term, "C42B Human ProstateCancer", "C42B_(Prostate_Cancer) Human")) %>% 
  mutate(Term = str_replace(Term, "H9 Human BoneMarrow Lymphoma", "H9_(Bone_Marrow_Lymphoma) Human")) %>% 
  mutate(Term = str_replace(Term, "WistarRat Hippocampus Stress", "Hippocampus_(Stress) WistarRat")) %>% 
  mutate(Term = str_replace(Term, "WistarRat Hippocampus", "Hippocampus WistarRat")) %>% 
  mutate(Term = str_replace(Term, "A549 Human Lung Carcinoma", "A549_(Lung_Carcinoma) Human")) %>% 
  mutate(Term = str_replace(Term, "MC3T3E1 Mouse Bone", "MC3T3E1_(Bone) Mouse")) %>% 
  mutate(TF = str_split_fixed(Term, " ", 5)[, 1]) %>% 
  mutate(PMID = str_split_fixed(Term, " ", 5)[, 2]) %>% 
  mutate(method = str_split_fixed(Term, " ", 5)[, 3]) %>% 
  mutate(tissue = str_split_fixed(Term, " ", 5)[, 4]) %>% 
  mutate(species = str_split_fixed(Term, " ", 5)[, 5]) %>% 
  group_by(cluster) %>%
  slice_head(n = 3) %>% 
  mutate(top = row_number()) %>% 
  select(c(cluster, Term, top, TF, PMID, tissue, species)) -> to_text_plot




to_text_plot %>% 
  as.data.frame() %>% 
  mutate(top = paste0(top, ".")) %>% 
  mutate(tissue = str_replace_all(tissue, "_", " ")) %>% 
  mutate(term = sprintf("%-2s %-7s %-29s %-2s", top, TF, tissue, PMID)) -> to_text_plot




# Przypisanie ID w taki sposób, by każde z nich odpowiadało jednemu z trzech wierszy w klastrze
to_text_plot <- to_text_plot %>%
  mutate(id = rep(c(1,2.5,4),18)) %>% 
  mutate(PMID = paste0("PMID:", PMID)) %>% 
  mutate(TF = paste(top, TF))
  


ggplot(to_text_plot, aes(x = 1, y = id, label = term)) +
  # geom_text(hjust = 0) +
  geom_text(aes(label = TF), hjust = 0) +  # Add term text aligned to the left
  geom_text(aes(label = tissue), hjust = 0, nudge_x = 0.06) +  # Add PMID text, aligned to the right
  geom_text(aes(label = species), hjust = 0, nudge_x = 0.21) +
  geom_text(aes(label = PMID), hjust = 0, nudge_x = 0.26) +
  facet_grid(cluster ~ .) +  # Use free_x to adjust x-axis per facet
  scale_x_continuous(limits = c(1, 1.34)) +  # Set y-axis limits from 0 to 4
  scale_y_reverse(limits = c(5, 0)) + 
  # scale_x_continuous(limits = c(1, 1.5), labels = c("Left", "Right")) + 
  theme(axis.text.x = element_blank(),  # Hide x-axis text
        axis.ticks.x = element_blank(),  # Hide x-axis ticks
        axis.title.x = element_blank(),  # Hide x-axis title
        axis.title.y = element_blank(),  # Hide y-axis title
        axis.text.y = element_blank(),   # Hide y-axis text
        axis.ticks.y = element_blank())  # Hide y-axis ticks


create_heatmap <- function(data, x_axis, y_axis, panel_spacing_size = 6, border_size = 1.5, tile_size = 0.7,
                           x_axis_text_size = 12, y_axis_text_size = 12) {
  # Sort y-axis elements in reverse alphabetical order and adjust as a factor
  data <- data %>%
    arrange(desc(!!sym(y_axis))) %>%
    mutate(!!y_axis := factor(!!sym(y_axis), levels = unique(!!sym(y_axis))))
  
  # Create the heatmap
  plot <- ggplot(data, aes(x = get(x_axis), y = get(y_axis), fill = factor(value))) +
    geom_tile(color = "black", size = border_size, width = tile_size, height = tile_size, na.rm = TRUE) + # Tiles with custom sizes
    scale_fill_manual(values = c("0" = "#e2e2e2", "1" = "#ffa590", "2" = "#ff4122", "3" = "#c61a09")) + # Custom colors
    # scale_fill_manual(values = c("0" = "grey", "1" = "orange", "2" = "red")) + # Custom colors
    theme_minimal() + # Minimal theme
    theme(
      axis.text.x.top = element_text(angle = 90, hjust = 1, size = x_axis_text_size, color = "black", vjust = 0.5), # X-axis labels at the top
      axis.ticks.x.top = element_line(color = "black"), # Ensure ticks are at top if needed
      axis.text.x = element_blank(), # Hide the default x-axis text at bottom
      axis.text.y = element_text(size = y_axis_text_size), # Y-axis text size
      axis.title.x.top = element_blank(), # Ensure no title at the top x-axis
      axis.title.x = element_blank(), # Hide the bottom x-axis title if it exists
      axis.title.y = element_blank(), # Hide y-axis title
      axis.ticks.length = unit(0, "points"), # Remove tick marks
      panel.spacing = unit(panel_spacing_size, "lines"), # Custom spacing between cells
      plot.background = element_rect(fill = "white", colour = NA), # Transparent plot background
      panel.grid.major = element_blank(), # Remove major grid lines
      panel.grid.minor = element_blank(), # Remove minor grid lines
      strip.background = element_rect(fill = "white", colour = "black"), # Adjust the background of facet labels
      panel.border = element_blank(), # Remove default panel border
      legend.position = "bottom", # Move legend to bottom
      legend.box = "horizontal", # Align legend items horizontally
      legend.title = element_blank() # Remove legend title if desired
    ) +
    scale_x_discrete(position = "top")  # Move the x-axis labels to the top
  
  # Add a vertical line on the left side only
  plot + geom_vline(xintercept = 0.5, color = "black", size = 1)
}


create_heatmap(data=to_plot2, x_axis = "variable", y_axis = "cluster",
               border_size = 1, tile_size = 0.7, x_axis_text_size = 14, y_axis_text_size = 14) -> p1






svg("results/google-drive/enrichr/gr-dependent-transcriptional-pattern/chea_TF_heatmap2-continous.svg", width = 6, height = 16)
create_heatmap(data=to_plot, x_axis = "variable", y_axis = "cluster",
               border_size = 1, tile_size = 0.7, x_axis_text_size = 14, y_axis_text_size = 14)

dev.off()

svg("results/google-drive/enrichr/gr-dependent-transcriptional-pattern/chea_text_plot.svg", width = 5, height = 13)
text_plot
dev.off()


create_heatmap_continous <- function(data, x_axis, y_axis, panel_spacing_size = 6, border_size = 1.5, tile_size = 0.7,
                           x_axis_text_size = 12, y_axis_text_size = 12, p_value_column) {
  # # Sort y-axis elements in reverse alphabetical order and adjust as a factor
  data <- data %>%
    arrange(desc(!!sym(y_axis))) %>%
    mutate(!!y_axis := factor(!!sym(y_axis), levels = unique(!!sym(y_axis))))
  print(data)
  # Create the heatmap
  plot <- ggplot(data, aes(x = get(x_axis), y = get(y_axis), fill = abs(value))) +
    geom_tile(color = "black", size = border_size, width = tile_size, height = tile_size, na.rm = TRUE) +
    scale_fill_gradient(low = "white",high = "#e68a89", na.value = "grey50") +# Tiles with custom sizes
  #   "white"             = "#f1bcbb",
  # "pastel_orange"     = "#edacab",
  # "pastel_red"        = "#e68a89"
    # scale_fill_gradient2(low = "white", mid = "red", high = "darkred") + # Continuous gradient based on p-value
    theme_minimal() + # Minimal theme
    theme(
      axis.text.x.top = element_text(angle = 90, hjust = 1, size = x_axis_text_size, color = "black", vjust = 0.5), # X-axis labels at the top
      axis.ticks.x.top = element_line(color = "black"), # Ensure ticks are at top if needed
      axis.text.x = element_blank(), # Hide the default x-axis text at bottom
      axis.text.y = element_text(size = y_axis_text_size), # Y-axis text size
      axis.title.x.top = element_blank(), # Ensure no title at the top x-axis
      axis.title.x = element_blank(), # Hide the bottom x-axis title if it exists
      axis.title.y = element_blank(), # Hide y-axis title
      axis.ticks.length = unit(0, "points"), # Remove tick marks
      panel.spacing = unit(panel_spacing_size, "lines"), # Custom spacing between cells
      plot.background = element_rect(fill = "white", colour = NA), # Transparent plot background
      panel.grid.major = element_blank(), # Remove major grid lines
      panel.grid.minor = element_blank(), # Remove minor grid lines
      strip.background = element_rect(fill = "white", colour = "black"), # Adjust the background of facet labels
      panel.border = element_blank(), # Remove default panel border
      legend.position = "bottom", # Move legend to bottom
      legend.box = "horizontal", # Align legend items horizontally
      legend.title = element_blank() # Remove legend title if desired
    )
    scale_x_discrete(position = "top")  # Move the x-axis labels to the top
  plot
  # Add a vertical line on the left side only
  plot + geom_vline(xintercept = 0.5, color = "black", size = 1)
}


create_heatmap_continous(data=to_plot3, x_axis = "variable", y_axis = "cluster",
               border_size = 1, tile_size = 0.7, x_axis_text_size = 14, y_axis_text_size = 14)

svg("results/google-drive/enrichr/gr-dependent-transcriptional-pattern/chea_TF_heatmap2-continous.svg", width = 6, height = 18)
create_heatmap(data=to_plot2, x_axis = "variable", y_axis = "cluster",
               border_size = 1, tile_size = 0.7, x_axis_text_size = 14, y_axis_text_size = 14)

dev.off()
