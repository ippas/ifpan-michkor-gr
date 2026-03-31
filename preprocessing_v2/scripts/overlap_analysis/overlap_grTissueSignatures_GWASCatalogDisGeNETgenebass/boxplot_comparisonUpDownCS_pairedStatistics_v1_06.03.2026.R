GWASCatalogDisGeNETgenebass_GrTissueSystemic_overlapChi2$processed$original_data$df %>% 
  mutate(
    log10_pvalue = -log10(p_value),
    log10_geneOverlap = log10(gene_overlap_count + 1),
    combine_score = log10_pvalue * log10_geneOverlap
  ) %>%
  filter(!(Var1 %in% c(
    "DisGeNET_NA", "genebass_NA", "GWASCatalog_NA",
    "genebass_F9x_Date_F99_first_reported_(mental_disorder,_not_otherwise_specified)"
  ))) %>%
  rowwise() %>%
  mutate(
    matched = list(str_extract_all(Var1, paste(patterns, collapse="|"))[[1]]),
    n_unique_patterns = n_distinct(matched)
  ) %>%
  ungroup() %>%
  filter(n_unique_patterns == 1) %>% 
  mutate(
    icd10_category = map_chr(matched, ~ .x[[1]])
  ) %>%
  # filter(!is.na(regulation), !is.na(icd10_category)) %>%
  filter(!is.na(icd10_category)) %>% 
  select(-c(log2_odds_ratio, fdr_value, fdr)) %>% 
  rename(
    grSignature        = Var2,
    phenotype          = Var1,
    observed_overlap   = gene_overlap_count,
    grSignature_nGenes = Var2_n_genes,
    phenotype_nGenes   = Var1_n_genes,
    mapped_icd10_range = icd10_category
  ) %>% 
  mutate(
    grSignature = case_when(
      str_detect(grSignature, "BloodCellsUp") ~ "bloodUp",
      str_detect(grSignature, "BloodCellsDown") ~ "bloodDown",
      str_detect(grSignature, "LungCellsUp") ~ "lungUp",
      str_detect(grSignature, "LungCellsDown") ~ "lungDown",
      str_detect(grSignature, "NeuralCellsUp") ~ "neuralUp",
      str_detect(grSignature, "NeuralCellsDown") ~ "neuralDown",
      str_detect(grSignature, "global_GR_genes_globalUp") ~ "systemicUp",
      str_detect(grSignature, "global_GR_genes_globalDown") ~ "systemicDown",
      TRUE ~ grSignature
    )
  ) %>% 
  mutate(
    # ---- regulation ----
    regulation = case_when(
      str_detect(grSignature, "Up") ~ "up",
      str_detect(grSignature, "Down") ~ "down",
      TRUE ~ NA_character_
    ),
    
    # ---- tissue ----
    tissue = case_when(
      str_detect(grSignature, "blood") ~ "blood",
      str_detect(grSignature, "lung") ~ "lung",
      str_detect(grSignature, "neural") ~ "neural",
      str_detect(grSignature, "systemic") ~ "systemic",
      TRUE ~ NA_character_
    ),
    
    # ---- signature type ----
    signatureType = case_when(
      tissue == "systemic" ~ "systemic",
      tissue %in% c("blood","lung","neural") ~ "tissue",
      TRUE ~ NA_character_
    )
  ) %>% 
  select(
    mapped_icd10_range,
    phenotype,
    phenotype_nGenes,
    grSignature,
    tissue,
    signatureType,
    regulation,
    grSignature_nGenes,
    expected_overlap,
    observed_overlap,
    odds_ratio,
    chi2,
    p_value,
    combine_score,
    overlap_genes,
    matched,
    n_unique_patterns
  ) %>% 
  select(-matched, -n_unique_patterns) -> df_preprocessing


df_preprocessing %>% 
  filter(p_value < 0.05 & observed_overlap >= 3) %>% 
  filter(tissue == "systemic") %>% .$phenotype %>% 
  unique -> phenotypes_vector


df_preprocessing %>% 
  filter(tissue == "systemic") %>% 
  filter(phenotype %in% phenotypes_vector) %>% 
  select(phenotype, regulation, combine_score, mapped_icd10_range, tissue, p_value, observed_overlap, overlap_genes) %>% 
  mutate(
    source    = str_extract(phenotype, "^[^_]+"),
    phenotype = str_remove(phenotype, "^[^_]+_F[0-9]x_")
  ) %>% 
  mutate(
    mapped_icd10_range = icd10_relabel(mapped_icd10_range)
  ) %>% 
  select(phenotype, source, mapped_icd10_range, regulation, tissue, combine_score, p_value, observed_overlap, overlap_genes) -> systemic_df

systemic_df %>%
  select(phenotype, regulation, combine_score) %>%
  pivot_wider(
    names_from = regulation,
    values_from = combine_score
  ) %>% 
  { wilcox.test(.$up, .$down, paired = TRUE) }


# ##############################################################################
df_preprocessing %>% 
  filter(p_value < 0.05 & observed_overlap >= 3) %>% 
  filter(tissue == "neural") %>% .$phenotype %>% 
  unique -> phenotypes_vector

df_preprocessing %>% 
  filter(tissue == "neural") %>% 
  filter(phenotype %in% phenotypes_vector) %>% 
  select(phenotype, regulation, combine_score, mapped_icd10_range, tissue, p_value, observed_overlap, overlap_genes) %>% 
  mutate(
    source    = str_extract(phenotype, "^[^_]+"),
    phenotype = str_remove(phenotype, "^[^_]+_F[0-9]x_")
  ) %>% 
  mutate(
    mapped_icd10_range = icd10_relabel(mapped_icd10_range)
  ) %>% 
  select(phenotype, source, mapped_icd10_range, regulation, tissue, combine_score, p_value, observed_overlap, overlap_genes) -> neural_df

neural_df %>%
  select(phenotype, regulation, combine_score) %>%
  pivot_wider(
    names_from = regulation,
    values_from = combine_score
  ) %>% 
  { wilcox.test(.$up, .$down, paired = TRUE) }

# ##############################################################################
df_preprocessing %>% 
  filter(p_value < 0.05 & observed_overlap >= 3) %>% 
  filter(tissue == "blood") %>% .$phenotype %>% 
  unique -> phenotypes_vector

df_preprocessing %>% 
  filter(tissue == "blood") %>% 
  filter(phenotype %in% phenotypes_vector) %>% 
  select(phenotype, regulation, combine_score, mapped_icd10_range, tissue, p_value, observed_overlap, overlap_genes) %>% 
  mutate(
    source    = str_extract(phenotype, "^[^_]+"),
    phenotype = str_remove(phenotype, "^[^_]+_F[0-9]x_")
  ) %>% 
  mutate(
    mapped_icd10_range = icd10_relabel(mapped_icd10_range)
  ) %>% 
  select(phenotype, source, mapped_icd10_range, regulation, tissue, combine_score, p_value, observed_overlap, overlap_genes) -> blood_df

blood_df %>%
  select(phenotype, regulation, combine_score) %>%
  pivot_wider(
    names_from = regulation,
    values_from = combine_score
  ) %>% 
  { wilcox.test(.$up, .$down, paired = TRUE) }
# ##############################################################################
df_preprocessing %>% 
  filter(p_value < 0.05 & observed_overlap >= 3) %>% 
  filter(tissue == "lung") %>% .$phenotype %>% 
  unique -> phenotypes_vector

df_preprocessing %>% 
  filter(tissue == "lung") %>% 
  filter(phenotype %in% phenotypes_vector) %>% 
  select(phenotype, regulation, combine_score, mapped_icd10_range, tissue, p_value, observed_overlap, overlap_genes) %>% 
  mutate(
    source    = str_extract(phenotype, "^[^_]+"),
    phenotype = str_remove(phenotype, "^[^_]+_F[0-9]x_")
  ) %>% 
  mutate(
    mapped_icd10_range = icd10_relabel(mapped_icd10_range)
  ) %>% 
  select(phenotype, source, mapped_icd10_range, regulation, tissue, combine_score, p_value, observed_overlap, overlap_genes) -> lung_df

lung_df %>%
  select(phenotype, regulation, combine_score) %>%
  pivot_wider(
    names_from = regulation,
    values_from = combine_score
  ) %>% 
  { wilcox.test(.$up, .$down, paired = TRUE) }

library(dplyr)
library(ggplot2)

rbind(
  systemic_df,
  neural_df,
  blood_df,
  lung_df
) %>%
  mutate(
    regulation = factor(regulation, levels = c("up", "down")),
    tissue    = factor(tissue, levels = c("systemic", "neural", "blood", "lung"))
  ) %>%
  # Najpewniejszy sposób – wymuszamy kolejność już na poziomie danych
  arrange(tissue, regulation) %>%   # ← dodane dla bezpieczeństwa
  
  ggplot(aes(x = regulation, y = combine_score)) +
  geom_boxplot(
    outlier.shape = NA,
    color = "black"
  ) +
  geom_jitter(
    aes(color = tissue),
    width = 0.25,
    alpha = 0.8,
    size = 2
  ) +
  facet_wrap(~ tissue, ncol = 4) +   # ncol = 4 → jeden rząd, kolejność z levels
  
  # Paski i gwiazdki – można uprościć (jeśli chcesz)
  geom_segment(
    data = . %>% distinct(tissue) %>% mutate(
      x = 1, xend = 1, y = 13.7, yend = 14
    ),
    aes(x = x, xend = xend, y = y, yend = yend),
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = . %>% distinct(tissue) %>% mutate(
      x = 1, xend = 2, y = 14, yend = 14
    ),
    aes(x = x, xend = xend, y = y, yend = yend),
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = . %>% distinct(tissue) %>% mutate(
      x = 2, xend = 2, y = 13.7, yend = 14
    ),
    aes(x = x, xend = xend, y = y, yend = yend),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = . %>% distinct(tissue) %>% mutate(
      x = 1.5,
      y = 14.15,
      label = c("*", "*", "n.s.", "*")   # kolejność: systemic, neural, blood, lung
    ),
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 6
  ) +
  
  scale_color_manual(values = c(
    systemic = "#335C67",
    neural   = "#E09F3E",
    blood    = "#9E2A2B",
    lung     = "#540B0E"
  )) +
  
  coord_cartesian(ylim = c(0, 15)) +
  
  labs(
    x = "regulation",
    y = "combine_score"
  ) +
  
  theme_classic() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 12, face = "bold")   # opcjonalnie – wyraźniejsze etykiety facetów
  ) -> p



svg(
  "/home/mateusz/projects/ifpan-michkor-gr/results_v2/overlap/overlap_grTissues_GWASCatalogDisGeNETgenebass/figures/boxplot/boxplot_wilcoxon_UpDownComparison_v1_08.03.2026.svg",
  width = 9,
  height = 5
)

print(p)

dev.off()


# ##############################################################################
# ---- prepare supplementary file ----
# ##############################################################################
# ---- output path ----
out_xlsx <- "/home/mateusz/projects/ifpan-michkor-gr/results_v2/overlap/overlap_grTissues_GWASCatalogDisGeNETgenebass/tables/supplementFile_pairedWilcoxonUpDown_v1_09.03.2026.xlsx"

# ---- sheet order ----
sheet_order <- c("systemic", "neural", "blood", "lung")

# ---- prepare data ----
supplement_df_list <- rbind(
  systemic_df,
  neural_df,
  blood_df,
  lung_df
) %>% 
  mutate(
    overlap_genes = str_replace_all(overlap_genes, ",", ", "),
    tissue = factor(tissue, levels = sheet_order)
  ) %>% 
  rename(
    number_overlap_genes = observed_overlap
  ) %>% 
  pivot_wider(
    id_cols = c(phenotype, source, mapped_icd10_range, tissue),
    names_from = regulation,
    values_from = c(combine_score, p_value, number_overlap_genes, overlap_genes),
    names_glue = "{regulation}_{.value}"
  ) %>% 
  arrange(tissue, phenotype, source, mapped_icd10_range) %>% 
  split(.$tissue)

# zachowaj kolejność arkuszy
supplement_df_list <- supplement_df_list[sheet_order]

# usuń ewentualne puste elementy
supplement_df_list <- supplement_df_list[!sapply(supplement_df_list, is.null)]

# ---- workbook ----
wb <- createWorkbook()

# ---- styles ----
header_style <- createStyle(
  textDecoration = "bold",
  halign = "center",
  valign = "center"
)

sci_3dec_style <- createStyle(numFmt = "0.000E+00")

# ---- write each tissue as separate sheet ----
for (sheet_name in names(supplement_df_list)) {
  
  df_sheet <- supplement_df_list[[sheet_name]] %>%
    mutate(across(where(is.factor), as.character)) %>%
    select(-tissue)
  
  addWorksheet(wb, sheetName = sheet_name)
  
  writeData(
    wb,
    sheet = sheet_name,
    x = df_sheet,
    withFilter = TRUE,
    headerStyle = header_style
  )
  
  # zablokuj pierwszy wiersz i pierwszą kolumnę
  freezePane(
    wb,
    sheet = sheet_name,
    firstActiveRow = 2,
    firstActiveCol = 2
  )
  
  # autoszerokość kolumn
  setColWidths(
    wb,
    sheet = sheet_name,
    cols = 1:ncol(df_sheet),
    widths = "auto"
  )
  
  # kolumny do formatowania naukowego
  sci_cols <- which(names(df_sheet) %in% c(
    "up_combine_score",
    "down_combine_score",
    "up_p_value",
    "down_p_value"
  ))
  
  if (length(sci_cols) > 0 && nrow(df_sheet) > 0) {
    addStyle(
      wb,
      sheet = sheet_name,
      style = sci_3dec_style,
      rows = 2:(nrow(df_sheet) + 1),
      cols = sci_cols,
      gridExpand = TRUE,
      stack = TRUE
    )
  }
}

# ---- save workbook ----
saveWorkbook(wb, out_xlsx, overwrite = TRUE)
