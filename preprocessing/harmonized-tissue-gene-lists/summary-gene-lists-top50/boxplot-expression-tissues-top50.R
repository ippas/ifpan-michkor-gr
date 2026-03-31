################################################################################
# 1. tissue-cell - top50 genes per list

# prepare data
gene_list_log2ratio_top_50 %>% 
  select(-data) %>% 
  unnest(top_50_genes) %>% 
  ungroup %>% 
  filter(simple_tissue %in% c("lung", "blood", "brain")) %>% 
  mutate(simple_tissue2 = case_when(
    cell %in% c("macrophages", "hMDM", "mBMDM", "Monocyte", "THP-1") ~ "myeloid_derived_cells_(monocytes_macrophages)",
    cell %in% c("Bcell", "Tcell", "NKcell", "NALM6", "REH-overexpression-GCR") ~ "lymphoid_derived_(T_B_NK_cells)",
    cell %in% c("ASM", "ASM-siCEBPO", "ASM-siCNTR") ~ "ASM",
    cell %in% c("A549", "H1944", "H1975", "H2122", "BEAS-2B", "HBE", "pHBECs") ~ "epithelial_like_cells",
    cell %in% c("mglia", "astrocyte") ~ "glia",
    TRUE ~ "toRemove"
  )) %>% 
  filter(simple_tissue2 != "toRemove") %>% 
  mutate(simple_tissue = paste(simple_tissue, simple_tissue2, sep = "_")) %>% 
  group_by(simple_tissue, regulation) %>% 
  nest() %>% 
  mutate(n_list = map(data, ~ .x %>% .$label_regulation %>% unique %>% length)) %>% 
  mutate(n_genes_tissue = map(data, ~ .x %>% .$hgnc_symbol %>% unique %>% length)) %>% 
  unnest(c(n_list, n_genes_tissue)) -> top50_per_celltype


top50_plot_data <- rbind(plot_data, top50_per_celltype) %>%
  unnest(data) %>%
  group_by(simple_tissue) %>%
  nest() %>%
  mutate(n_records = map_int(data, nrow)) %>%
  arrange(desc(n_records)) %>%
  mutate(simple_tissue = factor(simple_tissue, levels = simple_tissue)) %>%
  unnest(data) %>%
  ungroup() %>%
  mutate(label = str_replace_all(label_regulation, "_down", ""),
         label = str_replace_all(label, "_up", "")) %>%
  group_by(simple_tissue) %>%
  nest() %>%
  mutate(n_expreriments = map_int(data, ~ .x$label %>% unique() %>% length),
         n_lists = map_int(data, ~ .x$label_regulation %>% unique() %>% length),
         n_genes_up = map_int(data, ~ .x %>% filter(regulation == "up") %>% pull(hgnc_symbol) %>% unique() %>% length),
         n_genes_down = map_int(data, ~ .x %>% filter(regulation == "down") %>% pull(hgnc_symbol) %>% unique() %>% length),
         n_records_up = map_int(data, ~ .x %>% filter(regulation == "up") %>% nrow()),
         n_records_down = map_int(data, ~ .x %>% filter(regulation == "down") %>% nrow())) %>%
  unnest(data) %>%
  ungroup() %>%
  mutate(label_figure = paste(simple_tissue, n_expreriments, n_lists, n_genes_up, n_genes_down, n_records, n_records_down, sep = "\n")) %>% 
  mutate(simple_tissue = factor(simple_tissue, levels = c(
    "lung",
    "lung_epithelial_like_cells",
    "lung_ASM",
    "blood",
    "blood_myeloid_derived_cells_(monocytes_macrophages)",
    "blood_lymphoid_derived_(T_B_NK_cells)",
    "bone",
    "brain",
    "brain_glia",
    "kidney",
    "adrenal-gland",
    "embryos",
    "muscle",
    "adipose",
    "cartilage",
    "liver",
    "small-intestine",
    "spleen",
    "placenta"
  ))) %>% 
  mutate(label_figure = factor(label_figure, levels = unique(label_figure[order(match(simple_tissue, levels(simple_tissue)))])))


# visualization data; top50 genes
svg("results/figures/harmonized-tissue-gene-lists/summary-gene-lists-top50/boxplot-log2ratio-top50-per-list-v2.svg", width = 8, height = 10)
ggplot(top50_plot_data, aes(x = label_figure, y = log2ratio, color = regulation)) +
  geom_jitter(alpha = 0.2, size = 0.5) +
  geom_boxplot(aes(group = interaction(simple_tissue, regulation)), 
               outlier.shape = NA, color = "black", fill = NA, size = 0.5,
               position = "identity") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(hjust = 0, vjust = 0.5),
    axis.title.x = element_blank(),
    axis.ticks.x.top = element_line(),
    axis.text.x.top = element_text(size = rel(1))
  ) +
  scale_x_discrete(position = "top") +
  scale_color_manual(values = c("down" = "blue4", "up" = "firebrick"))
dev.off()



################################################################################
# 2. tissue-cell signatures

signature_plot_data <- rbind(plot_data, top50_per_celltype) %>%
  unnest(data) %>%
  group_by(simple_tissue, regulation, hgnc_symbol) %>%
  summarise(sum_rs = sum(rank_score, na.rm = TRUE), .groups = "drop") %>%
  group_by(simple_tissue, regulation) %>%
  arrange(desc(sum_rs)) %>%
  slice_head(n = 50) %>%
  left_join(
    rbind(plot_data, top50_per_celltype) %>% unnest(data),
    by = c("simple_tissue", "regulation", "hgnc_symbol")
  ) %>%
  distinct(simple_tissue, regulation, hgnc_symbol, label_regulation, log2ratio) %>%
  mutate(
    label = str_replace_all(label_regulation, "_down|_up", "")
  ) %>%
  group_by(simple_tissue) %>%
  mutate(
    n_expreriments = n_distinct(label),
    n_lists = n_distinct(label_regulation),
    n_genes_up = n_distinct(hgnc_symbol[regulation == "up"]),
    n_genes_down = n_distinct(hgnc_symbol[regulation == "down"]),
    n_records = n(),
    n_records_down = sum(regulation == "down"),
    label_figure = paste(simple_tissue, n_expreriments, n_lists, n_genes_up, n_genes_down, n_records, n_records_down, sep = "\n")
  ) %>%
  ungroup() %>%
  mutate(
    simple_tissue = factor(simple_tissue, levels = c(
      "lung", "lung_epithelial_like_cells", "lung_ASM",
      "blood", "blood_myeloid_derived_cells_(monocytes_macrophages)",
      "blood_lymphoid_derived_(T_B_NK_cells)", "bone", "brain", "brain_glia",
      "kidney", "adrenal-gland", "embryos", "muscle", "adipose",
      "cartilage", "liver", "small-intestine", "spleen", "placenta"
    )),
    label_figure = factor(label_figure, levels = unique(label_figure[order(match(simple_tissue, levels(simple_tissue)))]))
  )

svg("results/figures/harmonized-tissue-gene-lists/summary-gene-lists-top50/boxplot-log2ratio-signatures-tissue-cells-v2.svg", width = 8, height = 10)

ggplot(signature_plot_data, aes(x = label_figure, y = log2ratio, color = regulation)) +
  geom_jitter(alpha = 0.6, size = 0.5) +
  geom_boxplot(aes(group = interaction(simple_tissue, regulation)),
               outlier.shape = NA, color = "black", fill = NA, size = 0.5,
               position = "identity") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(hjust = 0, vjust = 0.5),
    axis.title.x = element_blank(),
    axis.ticks.x.top = element_line(),
    axis.text.x.top = element_text(size = rel(1))
  ) +
  scale_x_discrete(position = "top") +
  scale_color_manual(values = c("down" = "blue4", "up" = "firebrick"))

dev.off()



# rm(top50_per_celltype, top50_plot_data, signature_plot_data)



