library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

# ---- 1) Zdefiniuj sygnatury jako listę (nazwa -> wektor genów) ----
sig_list <- list(
  Neural_Up   = flat_allGrSignatures_31.10.2025[sig_names]$minusGlobalUpDown5TissuesDerivedCells_NeuralCellsUp,
  Neural_Down = flat_allGrSignatures_31.10.2025[sig_names]$minusGlobalUpDown5TissuesDerivedCells_NeuralCellsDown,
  Lung_Up     = flat_allGrSignatures_31.10.2025[sig_names]$minusGlobalUpDown5TissuesDerivedCells_LungCellsUp,
  Lung_Down   = flat_allGrSignatures_31.10.2025[sig_names]$minusGlobalUpDown5TissuesDerivedCells_LungCellsDown,
  Blood_Up    = flat_allGrSignatures_31.10.2025[sig_names]$minusGlobalUpDown5TissuesDerivedCells_BloodCellsUp,
  Blood_Down  = flat_allGrSignatures_31.10.2025[sig_names]$minusGlobalUpDown5TissuesDerivedCells_BloodCellsDown
)

# ---- 2) Zdefiniuj odpowiadające tabele GTEx jako listę ----
gtex_list <- list(
  Neural_Up   = nuralCellUp_gtex,
  Neural_Down = nuralCellDown_gtex,
  Lung_Up     = lungCellUp_gtex,
  Lung_Down   = lungCellDown_gtex,
  Blood_Up    = bloodCellUp_gtex,
  Blood_Down  = bloodCellDown_gtex
)

tissues_keep <- c("Brain", "Whole_Blood", "Lung")
thresholds <- c(1, 5, 10)

# ---- 3) Helper: liczy podsumowanie dla jednej sygnatury ----
summarize_one_signature <- function(sig_name, sig_genes, gtex_df,
                                    tissues_keep = tissues_keep,
                                    thresholds = thresholds) {
  
  sig_genes <- unique(sig_genes)
  n_original <- length(sig_genes)
  
  gtex_df2 <- gtex_df %>%
    filter(tissue %in% tissues_keep) %>%
    # upewnij się, że geneSymbol jest unikalny w obrębie tissue (na wszelki)
    distinct(geneSymbol, tissue, .keep_all = TRUE)
  
  # globalnie: ile genów w ogóle wykrytych w GTEx (w którymkolwiek tissue z tissues_keep)
  n_detected_global <- gtex_df2 %>%
    distinct(geneSymbol) %>%
    nrow()
  
  # per tissue: ile wykrytych + % powyżej progów
  per_tissue <- gtex_df2 %>%
    group_by(tissue) %>%
    summarise(
      signature = sig_name,
      n_genes_original = n_original,
      n_genes_detected_in_gtex_global = n_detected_global,
      n_genes_detected_in_gtex_tissue = n(),
      pct_detected_in_tissue = 100 * n() / n_original,
      
      pct_gt_1  = 100 * mean(median_max > thresholds[1], na.rm = TRUE),
      pct_gt_5  = 100 * mean(median_max > thresholds[2], na.rm = TRUE),
      pct_gt_10 = 100 * mean(median_max > thresholds[3], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    relocate(signature, tissue)
  
  # jeśli w jakimś tissue nie ma żadnego genu, to dodaj wiersz z zerami
  per_tissue %>%
    right_join(tibble(tissue = tissues_keep), by = "tissue") %>%
    mutate(
      signature = sig_name,
      n_genes_original = n_original,
      n_genes_detected_in_gtex_global = n_detected_global,
      n_genes_detected_in_gtex_tissue = replace_na(n_genes_detected_in_gtex_tissue, 0L),
      pct_detected_in_tissue = replace_na(pct_detected_in_tissue, 0),
      pct_gt_1  = replace_na(pct_gt_1, 0),
      pct_gt_5  = replace_na(pct_gt_5, 0),
      pct_gt_10 = replace_na(pct_gt_10, 0)
    ) %>%
    arrange(tissue)
}

# ---- 4) Budowa tabeli zbiorczej ----
signature_summary_tbl <- imap_dfr(sig_list, ~{
  sig_name  <- .y
  sig_genes <- .x
  gtex_df   <- gtex_list[[sig_name]]
  summarize_one_signature(sig_name, sig_genes, gtex_df,
                          tissues_keep = tissues_keep,
                          thresholds = thresholds)
})

signature_summary_tbl


lungCellDown_gtex %>% 
  filter(tissue %in% c("Brain", "Lung")) %>% 
  filter(query_gene_symbol %in% c(
    "CDC25A", "CDH6", "CREB5", "EML1",
    "KCTD15", "MFHAS1", "MTCL1", "MYO10",
    "NRG1", "RPS6KA5", "SAMD5",
    "SEMA3E", "SLC25A37", "TRAF3", "VCAN"
  )) %>% as.data.frame() %>% 
  select(-c(unit, n_subtissues, median_mean, median_median, median_sd, median_min)) %>% 
  arrange(median_max)



lungCellDown_gtex %>% 
  filter(tissue %in% c("Brain", "Lung")) %>% 
  filter(query_gene_symbol %in% c(
    "CDC25A", "CDH6", "CREB5", "EML1",
    "KCTD15", "MFHAS1", "MTCL1", "MYO10",
    "NRG1", "RPS6KA5", "SAMD5",
    "SEMA3E", "SLC25A37", "TRAF3", "VCAN"
  )) %>% 
  ggplot(aes(x = tissue, y = median_max)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.8) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  theme_classic() +
  labs(
    x = "Tissue",
    y = "Median max expression (TPM)",
    title = "Comparison of median_max expression between Lung and Brain"
  )


genes <- c(
  "CDC25A", "CDH6", "CREB5", "EML1",
  "KCTD15", "MFHAS1", "MTCL1", "MYO10",
  "NRG1", "RPS6KA5", "SAMD5",
  "SEMA3E", "SLC25A37", "TRAF3", "VCAN"
)

# przygotowanie danych + sparowanie
df_plot <- lungCellDown_gtex %>% 
  filter(tissue %in% c("Brain", "Lung")) %>% 
  filter(query_gene_symbol %in% genes) %>% 
  select(query_gene_symbol, tissue, median_max) %>% 
  mutate(log2_median_max = log2(median_max + 1))

df_paired <- df_plot %>%
  select(query_gene_symbol, tissue, log2_median_max) %>%
  pivot_wider(
    names_from = tissue,
    values_from = log2_median_max
  ) %>%
  drop_na(Brain, Lung)

# paired t-test
t_res <- t.test(df_paired$Brain, df_paired$Lung, paired = TRUE)

pval <- signif(t_res$p.value, 3)

# wykres
ggplot(df_plot, aes(x = tissue, y = log2_median_max)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.8) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  theme_classic() +
  labs(
    x = "Tissue",
    y = "log2(Median max expression + 1) [TPM]",
    title = "Paired comparison of median_max expression between Lung and Brain",
    subtitle = paste0("Paired t-test p-value = ", pval)
  )
