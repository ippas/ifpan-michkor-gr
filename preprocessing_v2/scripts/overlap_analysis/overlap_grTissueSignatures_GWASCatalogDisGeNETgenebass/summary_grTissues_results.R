grTissue_summaryResults <- data.frame(
  F_category = c("F0x","F1x","F2x","F3x","F4x","F5x","F6x","F7x","F8x","F9x"),
  
  n_association = c(5,9,5,23,6,1,0,0,2,2),
  n_genes       = c(20,53,62,114,49,10,0,0,13,16),
  
  brain_n_association = c(2,3,2,7,3,1,0,0,0,0),
  brain_n_genes       = c(6,9,27,28,28,10,0,0,0,0),
  
  blood_n_association = c(2,4,1,4,1,0,0,0,0,0),
  blood_n_genes       = c(10,23,10,20,9,0,0,0,0,0),
  
  lung_n_association = c(1,2,2,12,2,0,0,0,2,2),
  lung_n_genes       = c(4,21,26,68,12,0,0,0,13,16),
  
  all_phenotypes = c(24,44,27,70,27,12,4,11,17,8),
  all_genes      = c(4474,6157,5555,7558,4973,3295,1381,1539,2834,1099)
) %>% 
  mutate(
    norm_assoc = n_association / all_phenotypes,
    norm_genes = n_genes / all_genes,
    
    brain_norm_assoc = brain_n_association / all_phenotypes,
    brain_norm_genes = brain_n_genes / all_genes,
    
    blood_norm_assoc = blood_n_association / all_phenotypes,
    blood_norm_genes = blood_n_genes / all_genes,
    
    lung_norm_assoc = lung_n_association / all_phenotypes,
    lung_norm_genes = lung_n_genes / all_genes
  )

df_long <- grTissue_summaryResults %>%
  select(F_category, ends_with("norm_genes")) %>%
  pivot_longer(
    cols = -F_category,
    names_to = "tissue",
    values_to = "norm_genes"
  ) %>%
  mutate(
    tissue = recode(tissue,
                    norm_genes = "All",
                    brain_norm_genes = "Brain",
                    blood_norm_genes = "Blood",
                    lung_norm_genes  = "Lung"
    )
  )

# --- sortowanie wg ALL malejąco ---
order_levels <- df_long %>%
  filter(tissue == "All") %>%
  arrange(norm_genes) %>%
  pull(F_category)

df_long <- df_long %>%
  mutate(
    F_category = factor(F_category, levels = order_levels),
    tissue = factor(tissue, levels = c("All","Brain","Blood","Lung"))
  )


ggplot(df_long, aes(x = norm_genes, y = F_category)) +
  geom_col(width = 0.7, fill = "#8B5A2B", color = "black", size = 1) +
  facet_wrap(~tissue, ncol = 4) +
  theme_classic() +
  labs(
    x = "Fraction of associated genes",
    y = NULL
  )

GWASCatalogDisGeNETgenebass_GrTissue_overlapChi2$processed$original_data$df %>%
  filter(grepl("_F", Var1)) %>%
  filter(p_value < 0.05) %>%
  filter(gene_overlap_count > 2) %>%
  filter(log2_odds_ratio > 0) %>% 
  mutate(
    F_category = factor(
      str_extract(Var1, "F[0-9]x"),
      levels = paste0("F", 0:9, "x")
    )
  ) %>%
  filter(Var1 != "GWASCatalog_F9x_F3x_F8x_F2x_attention deficit hyperactivity disorder,bipolar disorder,autism spectrum disorder,schizophrenia,major depressive disorder") %>% 
  filter(Var1 != "genebass_F9x_Date_F99_first_reported_(mental_disorder,_not_otherwise_specified)") %>% 
  na.omit() -> GWASCatalogDisGeNETgenebass_GrTissue_overlapChi2_dfFilter


GWASCatalogDisGeNETgenebass_GrTissue_overlapChi2_dfFilter %>% filter(grepl("neural", Var2, ignore.case = T)) %>% 
  .$overlap_genes %>% strsplit(",") %>% 
  unlist %>% 
  unique()

multi_downloadMedianSubtissuesGTEx_tissueSummaryExpression(gene_symbol = GWASCatalogDisGeNETgenebass_GrTissue_overlapChi2_dfFilter %>% filter(grepl("neural", Var2, ignore.case = T)) %>% 
                                                             .$overlap_genes %>% strsplit(",") %>% 
                                                             unlist %>% 
                                                             unique()) -> brainGenesOverlap_gtex

multi_downloadMedianSubtissuesGTEx_tissueSummaryExpression(gene_symbol = GWASCatalogDisGeNETgenebass_GrTissue_overlapChi2_dfFilter %>% filter(grepl("lung", Var2, ignore.case = T)) %>% 
                                                             .$overlap_genes %>% strsplit(",") %>% 
                                                             unlist %>% 
                                                             unique()) -> lungGenesOverlap_gtex


multi_downloadMedianSubtissuesGTEx_tissueSummaryExpression(gene_symbol = GWASCatalogDisGeNETgenebass_GrTissue_overlapChi2_dfFilter %>% filter(grepl("blood", Var2, ignore.case = T)) %>% 
                                                             .$overlap_genes %>% strsplit(",") %>% 
                                                             unlist %>% 
                                                             unique()) -> bloodGenesOverlap_gtex


brainGenesOverlap_gtex %>% 
  select(geneSymbol, tissue, median_max) %>% 
  filter(tissue %in% c("Brain", "Whole_Blood", "Lung"))
brainGenesOverlap_gtex
lungGenesOverlap_gtex
bloodGenesOverlap_gtex %>% 
  select(geneSymbol, tissue, median_max) %>% 
  filter(tissue %in% c("Brain", "Whole_Blood", "Lung")) %>%
  mutate(tissue = factor(tissue, levels = c("Brain", "Lung", "Whole_Blood"))) %>%
  ggplot(aes(x = tissue, y = median_max, fill = tissue)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.25) +
  scale_y_log10() +
  theme_classic() +
  labs(x = NULL, y = "GTEx median expression") +
  theme(legend.position = "none")


plot_df <-
  bind_rows(
    brainGenesOverlap_gtex  %>% mutate(gene_set = "Brain genes"),
    lungGenesOverlap_gtex   %>% mutate(gene_set = "Lung genes"),
    bloodGenesOverlap_gtex  %>% mutate(gene_set = "Blood genes")
  ) %>%
  select(geneSymbol, tissue, median_max, gene_set) %>%
  filter(tissue %in% c("Brain", "Whole_Blood", "Lung")) %>%
  mutate(
    tissue = factor(tissue, levels = c("Brain", "Lung", "Whole_Blood")),
    gene_set = factor(gene_set, levels = c("Brain genes","Lung genes","Blood genes"))
  )

ggplot(plot_df, aes(x = tissue, y = median_max, fill = tissue)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.25) +
  scale_y_log10() +
  facet_wrap(~gene_set, nrow = 1) +
  theme_classic() +
  labs(x = NULL, y = "GTEx median expression") +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

