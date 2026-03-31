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
  select(-split_parts) %>% head
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


overlap_results_preprocessing %>% 
  mutate(fdr = map(data, ~ .x$original_data$df$fdr)) %>% 
  mutate(phenotypes = map(data, ~ .x$original_data$df$Var1)) %>%
  mutate(gr_list = map(data, ~ .x$original_data$df$Var2)) %>% 
  mutate(gene_overlap_count = map(data, ~ .x$original_data$df$gene_overlap_count)) %>% 
  mutate(n_genes_category = map(data, ~{genes_phenotypes_PanUkBiobank[.x$original_data$rows] %>% unname() %>% unlist() %>% unique %>% length()})) %>% 
  unnest(n_genes_category)
  

tmp2$original_data$rows

genes_phenotypes_PanUkBiobank[tmp2$original_data$rows] %>% unname %>% unlist %>% unique %>% length()

