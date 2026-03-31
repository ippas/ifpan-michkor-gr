overlap_results_preprocessing %>% 
  mutate(fdr = map(data, ~ .x$original_data$df$fdr)) %>% 
  mutate(phenotypes = map(data, ~ .x$original_data$df$Var1)) %>%
  mutate(gr_list = map(data, ~ .x$original_data$df$Var2)) %>% 
  mutate(gene_overlap_count = map(data, ~ .x$original_data$df$gene_overlap_count)) %>% 
  mutate(n_genes_category = map(data, ~{genes_phenotypes_PanUkBiobank[.x$original_data$rows] %>% unname() %>% unlist() %>% unique %>% length()})) %>% 
  select(-data) %>% 
  # select(category, fdr, FDR_threshold, gene_overlap_count) %>% head
  unnest(c(fdr, gene_overlap_count, phenotypes, gr_list)) %>% 
  mutate(FDR_threshold_colors = ifelse(FDR_threshold == "<0.001", "#F8766D",
                                       ifelse(FDR_threshold == ">0.05", "#619CFF", "#00BA38"))) %>%
  mutate(FDR_threshold_colors = ifelse(gene_overlap_count < 3, "gray", FDR_threshold_colors)) %>% 
  mutate(category = as.factor(category)) %>%
  mutate(log10_fdr = -log10(fdr)) %>% 
  # filter(log10_fdr >  2) %>% as.data.frame() %>% 
  filter(FDR_threshold_colors != "gray") %>% 
  mutate(category2 = paste(category, " (", number_phenotypes, ")", sep = "")) -> tmp


# Corrected color map with standard color names or hexadecimal codes
color_tissue_map <- c(
  "brain" = "#f9dcdc", # Standard color name
  "adrenal-gland" = "#FFFF00", # Hexadecimal for bright yellow
  "skeletal-tissue" = "ivory", # A close match to "off-white"
  "embryos" = "#D1FFBD",
  "blood" = "#ff8282",
  "lung" = "lightblue",
  "kidney" = "#cc3333",
  "adipocyte-tissue" = "orange",
  "liver" = "saddlebrown", # A close match to "dark brown"
  "other" = "slategrey", # Correct spelling for slate gray
  "gray" = "gray"
)

color_tissue_map <- c(
  brain = "#E9B7C4",  # Soft pink-gray
  `adrenal-gland` = "#FFA500",  # Vibrant orange
  `skeletal-tissue` = "#EEE8AA",  # Soft beige
  embryos = "#98FB98",  # Pastel green
  blood = "#8B0000",  # Deep red
  lung = "#ADD8E6",  # Light blue
  kidney = "#CD853F",  # Earthy brown-red
  `adipocyte-tissue` = "#FFF8DC",  # Creamy white
  liver = "#8B4513",  # Dark red-brown
  other = "#7B68EE"  # Neutral purple
)

tmp$gr_list %>%
  unique() %>%
  # Create a tibble with the original gr_list for later use
  tibble(gr_list = .) %>%
  # Split gr_list into parts and create new columns
  mutate(split_parts = strsplit(gr_list, split = "_")) %>%
  mutate(pmid = map_chr(split_parts, ~ .x[1]),
         tissue = map_chr(split_parts, ~ .x[2]),
         cell = map_chr(split_parts, ~ .x[3]),
         regulation = map_chr(split_parts, ~ .x[4])) %>%
  mutate(regulation = ifelse(gr_list == "pmid:24777604_embryos_hypothalamic-region_NPSCs|NPSCs-KO-Cav1_down", "down", regulation)) %>% 
  mutate(regulation = ifelse(gr_list == "pmid:24777604_embryos_hypothalamic-region_NPSCs|NPSCs-KO-Cav1_up", "up", regulation)) %>% 
  mutate(regulation = ifelse(gr_list == "pmid:26606517_embryos_hypothalamic-region_NPSCs_up", "up", regulation)) %>% 
  mutate(regulation = ifelse(gr_list == "pmid:26606517_embryos_hypothalamic-region_NPSCs_down", "down", regulation)) %>% 
  # Now you can safely remove split_parts column if it's no longer needed
  select(-split_parts) %>%
  mutate(tissue = ifelse(tissue == "NA", "brain", tissue)) %>%
  mutate(simple_tissue = case_when(
    tissue %in% c(
      "cortex", "hippocampus", "striatum", "prefrontal-cortex", "brain", 
      "hippocampal-progenitor-cell-line", "hippocampal-slices", "neuronal-pc12-cells", 
      "hippocampal-progenitor-cell-lie", "ARC", "hypothalamus", "pituitary-gland", 
      "cells-derived-from-human-fetal-brain-tissue"
    ) ~ "brain",
    tissue %in% c("kidney", "kidneys") ~ "kidney",
    tissue %in% c("adipose", "perigonadal-adipose-tissue", "adipocyte-tissue") ~ "adipocyte-tissue",
    tissue %in% c("blood", "blood-COPD") ~ "blood",
    tissue %in% c("skeletal-muscle", "anterior-thigh") ~ "skeletal-tissue",
    tissue %in% c("adrenal-cortex", "adrenal-gland") ~ "adrenal-gland",
    TRUE ~ tissue
  )) %>%
  group_by(simple_tissue) %>%
  mutate(n_simple_tissue = n()) %>% 
  ungroup() %>%
  mutate(simple_tissue2 = ifelse(n_simple_tissue <= 2, "other", simple_tissue)) %>%
  mutate(tissue_color = color_tissue_map[simple_tissue2]) %>% 
  mutate(regulation_color = case_when(
    regulation == "up" ~  "#FF5733",
    regulation == "down" ~ "#3498DB",
    TRUE ~ "#F1C40F" # Default case
  )) -> gr_list_metadata



left_join(tmp, gr_list_metadata, by = "gr_list") %>% as.data.frame() %>% 
  # filter(simple_tissue2 %in% c("brain", "adrenal-gland")) %>% 
  mutate(tissue_color = ifelse(fdr > 0.01, "gray", tissue_color),
         category2 = as.factor(category2),  # Assuming 'category2' needs to be factored
         log10_fdr = -log10(fdr)) %>%
  ggplot(aes(x = category2, y = log10_fdr, color = as.numeric(n_genes_category))) +
  # scale_color_identity(guide = "legend") +  # Use the actual colors specified in 'tissue_color'
  geom_jitter(alpha = 1) +
  # scale_color_manual(values = color_tissue_map) +
  labs(title = "Gene Overlap Manhattan Plot", x = "Category", y = "-log10(FDR)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  map(c(2, 10), ~geom_hline(yintercept = .x, linetype = "dashed", color = "black", size = 1))  +
  scale_color_manual(values = c("#E9B7C4" = "#E9B7C4",
                                "#FFA500" = "#FFA500",
                                "#EEE8AA" = "#EEE8AA",
                                "#98FB98" = "#98FB98",
                                "#8B0000" = "#8B0000",
                                "#ADD8E6" = "#ADD8E6",
                                "#CD853F" = "#CD853F",
                                "#FFF8DC" = "#FFF8DC",
                                "#8B4513" =  "#8B4513",
                                "#7B68EE" = "#7B68EE",
                                "gray" = "gray"),
                     labels = c("#E9B7C4" = "brain",
                                "#FFA500" = "Adrenal Gland",
                                "#EEE8AA" = "skeletal-tissue",
                                "#98FB98" = "embryos",
                                "#8B0000" = "blood",
                                "#ADD8E6" = "lung",
                                "#CD853F" = "kidney",
                                "#FFF8DC" = "adipocyte-tissue",
                                "#8B4513" = "liver",
                                "#7B68EE" = "other tissues",
                                "gray" = "n.s."),
                     guide = guide_legend(title = "Tissue Color"))
  

left_join(tmp, gr_list_metadata, by = "gr_list") %>% as.data.frame() %>% 
  # filter(simple_tissue2 %in% c("brain", "adrenal-gland")) %>% 
  mutate(tissue_color = ifelse(fdr > 0.01, "gray", tissue_color),
         category2 = as.factor(category2),  # Assuming 'category2' needs to be factored
         log10_fdr = -log10(fdr)) %>% 
  filter(fdr < 0.01) %>% 
  # filter(simple_tissue2 == "brain") %>% 
  select(c(category, simple_tissue2)) %>% 
  group_by(category) %>% 
  nest %>% 
  mutate(n_brain = map(data,  ~{.x$simple_tissue2  %>% table %>% as.data.frame() %>% set_colnames(c("tissue", "frequency"))})) %>% 
  # head() %>% 
  select(-data) %>% 
  unnest(n_brain) %>% filter(frequency > 2) %>% 
  filter(tissue == "skeletal-tissue") %>% as.data.frame()



left_join(tmp, gr_list_metadata, by = "gr_list") %>% as.data.frame() %>% 
  # filter(simple_tissue2 %in% c("brain", "adrenal-gland")) %>% 
  mutate(tissue_color = ifelse(fdr > 0.01, "gray", tissue_color),
         category2 = as.factor(category2),  # Assuming 'category2' needs to be factored
         log10_fdr = -log10(fdr)) %>% 
  filter(fdr < 0.01) %>% 
  filter(simple_tissue2 == "brain") %>% filter(fdr < 0.01) %>%  .$gr_list %>% unique()

left_join(tmp, gr_list_metadata, by = "gr_list") %>% 
  filter(fdr > 0.01) %>% as.data.frame() %>%
  filter(gene_overlap_count > 10) %>%
  # filter(gene_overlap_count < 10) %>% 
  ggplot(aes(x = category2, y = gene_overlap_count,  color = tissue_color)) +
  # geom_boxplot() +
  geom_jitter() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right") +
  # scale_color_identity()
  scale_color_manual(values = c("#E9B7C4" = "#E9B7C4",
                                "#FFA500" = "#FFA500",
                                "#EEE8AA" = "#EEE8AA",
                                "#98FB98" = "#98FB98",
                                "#8B0000" = "#8B0000",
                                "#ADD8E6" = "#ADD8E6",
                                "#CD853F" = "#CD853F",
                                "#FFF8DC" = "#FFF8DC",
                                "#8B4513" =  "#8B4513",
                                "#7B68EE" = "#7B68EE",
                                "gray" = "gray"),
                     labels = c("#E9B7C4" = "brain",
                                "#FFA500" = "Adrenal Gland",
                                "#EEE8AA" = "skeletal-tissue",
                                "#98FB98" = "embryos",
                                "#8B0000" = "blood",
                                "#ADD8E6" = "lung",
                                "#CD853F" = "kidney",
                                "#FFF8DC" = "adipocyte-tissue",
                                "#8B4513" = "liver",
                                "#7B68EE" = "other tissues",
                                "gray" = "n.s."),
                     guide = guide_legend(title = "Tissue Color")) 
################################################################################


# Main plot data preparation
plot_data <- tmp %>%
  left_join(gr_list_metadata, by = "gr_list") %>%
  as.data.frame() %>%
  mutate(tissue_color = ifelse(fdr > 0.01, "gray", tissue_color),
         category2 = as.factor(category2),
         log10_fdr = -log10(fdr))

ggplot(plot_data, aes(x = category2, y = log10_fdr, color = tissue_color)) +
  geom_jitter(alpha = 1) +
  geom_hline(yintercept = c(2, 10), linetype = "dashed", color = "black", size = 1) +
  scale_color_manual(values = c("#E9B7C4" = "#E9B7C4",
                                "#FFA500" = "#FFA500",
                                "#EEE8AA" = "#EEE8AA",
                                "#98FB98" = "#98FB98",
                                "#8B0000" = "#8B0000",
                                "#ADD8E6" = "#ADD8E6",
                                "#CD853F" = "#CD853F",
                                "#FFF8DC" = "#FFF8DC",
                                "#8B4513" =  "#8B4513",
                                "#7B68EE" = "#7B68EE",
                                "gray" = "gray"),
                     labels = c("brain", "Adrenal Gland", "skeletal-tissue", "embryos",
                                "blood", "lung", "kidney", "adipocyte-tissue", "liver", "other tissues", "n.s."),
                     guide = guide_legend(title = "Tissue Color")) +
  labs(title = "Gene Overlap Manhattan Plot", x = "Category", y = "-log10(FDR)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom") +
  geom_text(data = label_data, aes(label = "example", x = category2, y = y), vjust = 0, hjust = 0.5)

  color_tissue_map <- c(
    brain = "#E9B7C4",  # Soft pink-gray
    `adrenal-gland` = "#FFA500",  # Vibrant orange
    `skeletal-tissue` = "#EEE8AA",  # Soft beige
    embryos = "#98FB98",  # Pastel green
    blood = "#8B0000",  # Deep red
    lung = "#ADD8E6",  # Light blue
    kidney = "#CD853F",  # Earthy brown-red
    `adipocyte-tissue` = "#FFF8DC",  # Creamy white
    liver = "#8B4513",  # Dark red-brown
    other = "#7B68EE"  # Neutral purple
  )
tmp %>%
  left_join(gr_list_metadata, by = "gr_list") %>%
  # mutate(tissue_color = ifelse(fdr > 0.01, "other", as.character(tissue_color)),  # Use 'other' for non-significant
  #        tissue_color = factor(tissue_color, levels = names(color_tissue_map)),  # Ensure factor levels match map
  #        category2 = as.factor(category2),
  #        log10_fdr = -log10(fdr)) %>%
  ggplot(aes(x = category2, y = log10_fdr, color = tissue_color)) +
  geom_jitter(alpha = 1) +

  labs(title = "Gene Overlap Manhattan Plot", x = "Category", y = "-log10(FDR)", color = "Tissue Color") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom") +
  scale_color_manual(values = color_tissue_map)# Use custom color map






left_join(tmp, gr_list_metadata, by = "gr_list") %>% as.data.frame() %>% 
  mutate(tissue_color = ifelse(fdr > 0.01, "gray", tissue_color),
         category2 = as.factor(category2),  # Assuming 'category2' needs to be factored
         log10_fdr = -log10(fdr)) %>%
  ggplot(aes(x = category2, y = log10_fdr, color = tissue_color)) +
  geom_jitter(alpha = 1) +  # Assuming you want to keep the points slightly dispersed
  scale_color_identity() +  # Use the actual colors specified in 'tissue_color'
  labs(title = "Gene Overlap Manhattan Plot", x = "Category", y = "-log10(FDR)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom") +
  map(c(2, 10), ~geom_hline(yintercept = .x, linetype = "dashed", color = "black", size = 1)) +
  guides(color = guide_legend(title = "Tissue Color"))
  mutate(tissue_color = ifelse(fdr > 0.01, "gray", tissue_color)) %>%
  # mutate(regulation_color = ifelse(fdr > 0.01, "gray", regulation_color)) %>% 
  # filter(-log10(fdr) < 11) %>%
  gene_overlap_manhattan_plot(data = ., x = category2, color_column = "tissue_color")

gene_overlap_manhattan_plot <- function(data, x, color_column = NULL, title = "Gene Overlap Manhattan Plot", y_intercepts = c(2, 10), alpha = 1) {
  
  data <- data %>%
    mutate(category = as.factor(category),
           log10_fdr = -log10(fdr))
  
  # Initialize the ggplot with the basic setup
  p <- ggplot(data, aes(x = {{x}}, y = log10_fdr)) +
    geom_jitter(alpha = alpha) +
    theme_minimal() +
    labs(title = title, x = "Category", y = "-log10(FDR)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom") 
    purrr::map(y_intercepts, ~geom_hline(yintercept = .x, linetype = "dashed", color = "black", size = 1))
  
  p <- p + aes_string(color = color_column) + scale_color_identity() 
  # Add horizontal lines at specified y-intercepts
  p <- p + purrr::map(y_intercepts, ~geom_hline(yintercept = .x, linetype = "dashed", color = "black", size = 1))
  
  # Add text annotations for significant overlaps
  # p <- p + geom_text(aes(label = signif_overlap), vjust = -0.5, hjust = 0.5, check_overlap = TRUE)
  
  # Add the unique value annotation at the top of each category
  # p <- p + geom_text(data = signif_overlap, aes(x = {{x}}, y = max_log10_fdr + 0.5, label = signif_overlap), size = 3, vjust = 0, hjust = 0.5, check_overlap = TRUE)
  
  
  
  return(p)
}

