flat_allGrSignatures_31.10.2025[sig_names]

nuralCellUp_gtex <- multi_downloadMedianSubtissuesGTEx_tissueSummaryExpression(gene_symbol = flat_allGrSignatures_31.10.2025[sig_names]$minusGlobalUpDown5TissuesDerivedCells_NeuralCellsUp)
nuralCellDown_gtex <- multi_downloadMedianSubtissuesGTEx_tissueSummaryExpression(gene_symbol = flat_allGrSignatures_31.10.2025[sig_names]$minusGlobalUpDown5TissuesDerivedCells_NeuralCellsDown)

lungCellUp_gtex <- multi_downloadMedianSubtissuesGTEx_tissueSummaryExpression(gene_symbol = flat_allGrSignatures_31.10.2025[sig_names]$minusGlobalUpDown5TissuesDerivedCells_LungCellsUp)
lungCellDown_gtex <- multi_downloadMedianSubtissuesGTEx_tissueSummaryExpression(gene_symbol = flat_allGrSignatures_31.10.2025[sig_names]$minusGlobalUpDown5TissuesDerivedCells_LungCellsDown)

bloodCellUp_gtex <- multi_downloadMedianSubtissuesGTEx_tissueSummaryExpression(gene_symbol = flat_allGrSignatures_31.10.2025[sig_names]$minusGlobalUpDown5TissuesDerivedCells_BloodCellsUp)
bloodCellDown_gtex <- multi_downloadMedianSubtissuesGTEx_tissueSummaryExpression(gene_symbol = flat_allGrSignatures_31.10.2025[sig_names]$minusGlobalUpDown5TissuesDerivedCells_BloodCellsDown)


nuralCellUp_gtex
nuralCellDown_gtex


nuralCellDown_gtex %>% 
  filter(tissue %in% c("Brain", "Whole_Blood", "Lung")) %>% 
  select(geneSymbol, tissue, n_subtissues, median_max)


nuralCellUp_gtex %>% 
  filter(tissue %in% c("Brain", "Whole_Blood", "Lung")) %>% 
  ggplot(aes(x = tissue, y = log2(median_max))) +
  geom_boxplot(outlier.shape = 21, outlier.size = 2) +
  theme_minimal() +
  labs(
    x = "Tissue",
    y = "Median max expression",
    title = "Distribution of median_max across tissues"
  )

lungCellUp_gtex %>% 
  filter(tissue %in% c("Brain", "Whole_Blood", "Lung")) %>% 
  ggplot(aes(x = tissue, y = log2(median_max))) +
  geom_boxplot(outlier.shape = 21, outlier.size = 2) +
  theme_minimal() +
  labs(
    x = "Tissue",
    y = "Median max expression",
    title = "Distribution of median_max across tissues"
  )

bloodCellUp_gtex %>% 
  filter(tissue %in% c("Brain", "Whole_Blood", "Lung")) %>% 
  ggplot(aes(x = tissue, y = log2(median_max))) +
  geom_boxplot(outlier.shape = 21, outlier.size = 2) +
  theme_minimal() +
  labs(
    x = "Tissue",
    y = "Median max expression",
    title = "Distribution of median_max across tissues"
  )



nuralCellDown_gtex %>% 
  filter(tissue %in% c("Brain", "Whole_Blood", "Lung")) %>% 
  ggplot(aes(x = tissue, y = log2(median_max))) +
  geom_boxplot(outlier.shape = 21, outlier.size = 2) +
  theme_minimal() +
  labs(
    x = "Tissue",
    y = "Median max expression",
    title = "Distribution of median_max across tissues"
  )

lungCellDown_gtex %>% 
  filter(tissue %in% c("Brain", "Whole_Blood", "Lung")) %>% 
  ggplot(aes(x = tissue, y = log2(median_max))) +
  geom_boxplot(outlier.shape = 21, outlier.size = 2) +
  theme_minimal() +
  labs(
    x = "Tissue",
    y = "Median max expression",
    title = "Distribution of median_max across tissues"
  )

bloodCellDown_gtex %>% 
  filter(tissue %in% c("Brain", "Whole_Blood", "Lung")) %>% 
  ggplot(aes(x = tissue, y = log2(median_max))) +
  geom_boxplot(outlier.shape = 21, outlier.size = 2) +
  theme_minimal() +
  labs(
    x = "Tissue",
    y = "Median max expression",
    title = "Distribution of median_max across tissues"
  )


tissues_keep <- c("Brain", "Whole_Blood", "Lung")

combined_gtex <- bind_rows(
  nuralCellUp_gtex   %>% mutate(signature_set = "Neural", direction = "Up"),
  nuralCellDown_gtex %>% mutate(signature_set = "Neural", direction = "Down"),
  lungCellUp_gtex    %>% mutate(signature_set = "Lung",   direction = "Up"),
  lungCellDown_gtex  %>% mutate(signature_set = "Lung",   direction = "Down"),
  bloodCellUp_gtex   %>% mutate(signature_set = "Blood",  direction = "Up"),
  bloodCellDown_gtex %>% mutate(signature_set = "Blood",  direction = "Down")
) %>%
  filter(tissue %in% tissues_keep) %>%
  mutate(
    tissue        = factor(tissue, levels = tissues_keep),
    signature_set = factor(signature_set, levels = c("Neural", "Lung", "Blood")),
    direction     = factor(direction, levels = c("Up", "Down")),
    log2_median_max = log2(median_max)
  ) %>%
  select(geneSymbol, tissue, n_subtissues, median_max, log2_median_max, signature_set, direction)

ggplot(combined_gtex, aes(x = tissue, y = log2_median_max, fill = signature_set)) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.7,
    outlier.shape = 21,
    outlier.size = 1.6
  ) +
  facet_wrap(~ direction, nrow = 1) +
  theme_minimal() +
  labs(
    x = "Tissue",
    y = "log2(median_max)",
    fill = "Signature set",
    title = "GTEx expression across tissues (split by Up vs Down genes)"
  )


target_map <- tibble(
  signature_set = c("Neural", "Lung", "Blood"),
  target_tissue = c("Brain", "Lung", "Whole_Blood")
)

delta_tbl <- combined_gtex %>%
  left_join(target_map, by = "signature_set") %>%
  mutate(group = ifelse(tissue == target_tissue, "target", "non_target")) %>%
  group_by(direction, signature_set, group) %>%
  summarise(med = median(expr, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = group, values_from = med) %>%
  mutate(delta = target - non_target)

delta_tbl %>%
  ggplot(aes(x = signature_set, y = delta)) +
  geom_col(width = 0.7) +
  facet_wrap(~ direction) +
  labs(x = NULL, y = "Δ median(expression): target − non-target") +
  theme_classic()
