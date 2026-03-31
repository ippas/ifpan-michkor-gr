# ============================================================
# 🧬 Preparation and cleanup of GR-dependent gene signatures
# ============================================================

library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(openxlsx)
library(readr)

# ============================================================
# 1️⃣ Define primary tissue-specific signatures (UP/DOWN)
# ============================================================

tissue_up_signatures_genes <- list(
  NeuralCellsUp = gene_summaryList_brain$n10_short_time$summary_up$enrichr_summary_source$genes_3pub,
  BloodCellsUp  = gene_summaryList_blood$n10_short_time$summary_up$enrichr_summary_source$genes_3pub,
  LungCellsUp   = gene_summaryList_lung$n10_short_time$summary_up$enrichr_summary_source$genes_3pub
)

tissue_down_signatures_genes <- list(
  NeuralCellsDown = gene_summaryList_brain$n10_short_time$summary_down$enrichr_summary_source$genes_3pub,
  BloodCellsDown  = gene_summaryList_blood$n10_short_time$summary_down$enrichr_summary_source$genes_3pub,
  LungCellsDown   = gene_summaryList_lung$n10_short_time$summary_down$enrichr_summary_source$genes_3pub
)

# ============================================================
# 2️⃣ Identify repeated genes across tissues (UP/DOWN)
# ============================================================

repeated_up <- list(
  tissue_up_signatures_genes$NeuralCellsUp,
  tissue_up_signatures_genes$BloodCellsUp,
  tissue_up_signatures_genes$LungCellsUp
) %>%
  unlist() %>% table() %>% as.data.frame() %>%
  set_colnames(c("hgnc_symbol", "freq")) %>% filter(freq > 1) %>% pull(hgnc_symbol)

repeated_down <- list(
  tissue_down_signatures_genes$NeuralCellsDown,
  tissue_down_signatures_genes$BloodCellsDown,
  tissue_down_signatures_genes$LungCellsDown
) %>%
  unlist() %>% table() %>% as.data.frame() %>%
  set_colnames(c("hgnc_symbol", "freq")) %>% filter(freq > 1) %>% pull(hgnc_symbol)

repeatedTissueCells <- list(
  repeatedTissueCellsUp = repeated_up,
  repeatedTissueCellsDown = repeated_down
)

# ============================================================
# 3️⃣ Remove repeated genes from tissue signatures
# ============================================================

remove_repeated <- function(sig_list, repeated_vec) {
  lapply(sig_list, \(genes) genes[!genes %in% repeated_vec])
}

tissueCellsUpSignaturesMinusRepeated_genes   <- remove_repeated(tissue_up_signatures_genes, repeated_up)
tissueCellsDownSignaturesMinusRepeated_genes <- remove_repeated(tissue_down_signatures_genes, repeated_down)

# ============================================================
# 4️⃣ Remove genes belonging to specific clusters
# ============================================================

remove_cluster_genes <- function(sig_list, cluster_genes) {
  lapply(sig_list, \(genes) genes[!genes %in% cluster_genes])
}

signaturesMinusClusters_up <- list(
  minusClusterP  = remove_cluster_genes(tissue_up_signatures_genes, marpiech_clusters$cluster_P),
  minusClusterO  = remove_cluster_genes(tissue_up_signatures_genes, marpiech_clusters$cluster_O),
  minusClusterK  = remove_cluster_genes(tissue_up_signatures_genes, marpiech_clusters$cluster_K),
  minusClustersKPO = remove_cluster_genes(
    tissue_up_signatures_genes,
    unlist(marpiech_clusters[c("cluster_K", "cluster_P", "cluster_O")], use.names = FALSE)
  )
)

signaturesMinusClusters_down <- list(
  minusClusterD = remove_cluster_genes(
    tissue_down_signatures_genes,
    marpiech_clusters$cluster_D
  )
)

# ============================================================
# 5️⃣ Identify global GR-dependent genes (from full database)
# ============================================================

filtered_data <- papers_data_preprocessing %>%
  filter(
    !time %in% c(
      "240h","720h","240h_vs_720h","9weeks","3months","10days","14wekks","8weeks",
      "168h","240-528h","14weeks","3weeks","7weeks","672h","main_effect_of_treatments"
    ),
    !simple_tissue %in% c("Other", "placenta"),
    n_genes >= 10
  ) %>%
  mutate(simple_tissue_regulation = paste0(simple_tissue, "_", regulation))

get_genes_by_regulation_and_threshold <- function(df, regulation_type, threshold) {
  df %>%
    filter(regulation == regulation_type) %>%
    group_by(simple_tissue_regulation) %>%
    nest() %>%
    mutate(genes = map(data, \(x) unique(x$hgnc_symbol))) %>%
    pull(genes) %>%
    unlist() %>%
    table() %>%
    as.data.frame() %>%
    setNames(c("hgnc_symbol", "freq")) %>%
    filter(freq >= threshold) %>%
    arrange(desc(freq)) %>%
    pull(hgnc_symbol)
}

global_GR_genes <- list(
  globalUp4TissuesDerivedCells   = get_genes_by_regulation_and_threshold(filtered_data, "up",   4),
  globalUp5TissuesDerivedCells   = get_genes_by_regulation_and_threshold(filtered_data, "up",   5),
  globalUp6TissuesDerivedCells   = get_genes_by_regulation_and_threshold(filtered_data, "up",   6),
  globalDown4TissuesDerivedCells = get_genes_by_regulation_and_threshold(filtered_data, "down", 4),
  globalDown5TissuesDerivedCells = get_genes_by_regulation_and_threshold(filtered_data, "down", 5),
  globalDown6TissuesDerivedCells = get_genes_by_regulation_and_threshold(filtered_data, "down", 6)
)

# ============================================================
# 6️⃣ Remove global genes from tissue-specific signatures
# ============================================================

remove_global_genes <- function(sig_list, global_vec) {
  lapply(sig_list, \(genes) genes[!genes %in% global_vec])
}

signaturesMinusGlobal <- list(
  minusGlobalUp4TissuesDerivedCells   = remove_global_genes(tissue_up_signatures_genes,   global_GR_genes$globalUp4TissuesDerivedCells),
  minusGlobalDown4TissuesDerivedCells = remove_global_genes(tissue_down_signatures_genes, global_GR_genes$globalDown4TissuesDerivedCells),
  minusGlobalUp5TissuesDerivedCells   = remove_global_genes(tissue_up_signatures_genes,   global_GR_genes$globalUp5TissuesDerivedCells),
  minusGlobalDown5TissuesDerivedCells = remove_global_genes(tissue_down_signatures_genes, global_GR_genes$globalDown5TissuesDerivedCells),
  minusGlobalUp6TissuesDerivedCells   = remove_global_genes(tissue_up_signatures_genes,   global_GR_genes$globalUp6TissuesDerivedCells),
  minusGlobalDown6TissuesDerivedCells = remove_global_genes(tissue_down_signatures_genes, global_GR_genes$globalDown6TissuesDerivedCells)
)

# ============================================================
# 6️⃣➕ NEW: Remove both global UP & DOWN (5-tissue) from all tissue signatures
# ============================================================

global_updown_5_vec <- unique(c(
  global_GR_genes$globalUp5TissuesDerivedCells,
  global_GR_genes$globalDown5TissuesDerivedCells
))

signaturesMinusGlobalUpDown5TissuesDerivedCells <- list(
  minusGlobalUpDown5TissuesDerivedCells = c(
    remove_global_genes(tissue_up_signatures_genes,   global_updown_5_vec),
    remove_global_genes(tissue_down_signatures_genes, global_updown_5_vec)
  )
)

# ============================================================
# 7️⃣ Combine all objects into one structured list
# ============================================================

allGrSignatures_31.10.2025 <- list(
  tissue_up_signatures_genes = tissue_up_signatures_genes,
  tissue_down_signatures_genes = tissue_down_signatures_genes,
  tissueCellsUpSignaturesMinusRepeated_genes = tissueCellsUpSignaturesMinusRepeated_genes,
  tissueCellsDownSignaturesMinusRepeated_genes = tissueCellsDownSignaturesMinusRepeated_genes,
  repeatedTissueCells = repeatedTissueCells,
  signaturesMinusClusters_up = signaturesMinusClusters_up,
  signaturesMinusClusters_down = signaturesMinusClusters_down,
  global_GR_genes = global_GR_genes,
  signaturesMinusGlobal = signaturesMinusGlobal,
  signaturesMinusGlobalUpDown5TissuesDerivedCells = signaturesMinusGlobalUpDown5TissuesDerivedCells
)

# ============================================================
# 8️⃣ Flatten nested structure into a data frame
# ============================================================

flatten_named_list <- function(x, parent = character(), sep = ".") {
  if (is.list(x)) {
    purrr::imap(x, ~ flatten_named_list(.x, c(parent, .y), sep)) |> unlist(recursive = FALSE)
  } else if (is.atomic(x)) {
    nm <- paste(parent, collapse = sep)
    setNames(list(as.character(x)), nm)
  } else list()
}

flat_allGrSignatures_31.10.2025 <- flatten_named_list(allGrSignatures_31.10.2025)

names(flat_allGrSignatures_31.10.2025) <- sapply(names(flat_allGrSignatures_31.10.2025), function(nm) {
  parts <- unlist(strsplit(nm, "\\."))
  new_name <- if (length(parts) >= 2) {
    paste0(parts[length(parts) - 1], "_", parts[length(parts)])
  } else nm
  
  new_name <- gsub("^tissue_up_signatures_genes_", "", new_name)
  new_name <- gsub("^tissue_down_signatures_genes_", "", new_name)
  new_name
})

df_allGrSignatures_31.10.2025 <- map2_df(
  flat_allGrSignatures_31.10.2025,
  names(flat_allGrSignatures_31.10.2025),
  \(x, nm) tibble(hgnc_symbol = x, signature_name = nm)
)

# ============================================================
# 💾 Save all results
# ============================================================

save_path <- "results_v2/GR_signatures/GR_signatures_31.10.2025"
dir.create(save_path, recursive = TRUE, showWarnings = FALSE)

saveRDS(flat_allGrSignatures_31.10.2025, file = file.path(save_path, "flat_allGrSignatures_31.10.2025.rds"))
saveRDS(allGrSignatures_31.10.2025,      file = file.path(save_path, "allGrSignatures_31.10.2025.rds"))
saveRDS(df_allGrSignatures_31.10.2025,   file = file.path(save_path, "df_allGrSignatures_31.10.2025.rds"))

write_tsv(df_allGrSignatures_31.10.2025, file.path(save_path, "df_allGrSignatures_31.10.2025.tsv"))

wb <- createWorkbook()
addWorksheet(wb, "GR_signatures")
writeData(wb, "GR_signatures", df_allGrSignatures_31.10.2025)
saveWorkbook(wb, file.path(save_path, "df_allGrSignatures_31.10.2025.xlsx"), overwrite = TRUE)

cat("
✅ Saved GR signature objects to: results_v2/GR_signatures/GR_signatures_31.10.2025/
   ├── flat_allGrSignatures_31.10.2025.rds
   ├── allGrSignatures_31.10.2025.rds
   ├── df_allGrSignatures_31.10.2025.rds
   ├── df_allGrSignatures_31.10.2025.tsv
   └── df_allGrSignatures_31.10.2025.xlsx
")

# ============================================================
# 📘 Export UP/DOWN signatures to XLSX
# ============================================================

flat_sigs_up   <- flat_allGrSignatures_31.10.2025[grepl("Up",   names(flat_allGrSignatures_31.10.2025))]
flat_sigs_down <- flat_allGrSignatures_31.10.2025[grepl("Down", names(flat_allGrSignatures_31.10.2025))]

list_to_wide_df <- function(lst) {
  if (length(lst) == 0) return(data.frame())
  max_len <- max(lengths(lst))
  data.frame(lapply(lst, function(x) c(x, rep(NA, max_len - length(x)))), check.names = FALSE)
}

signaturesUp_df   <- list_to_wide_df(flat_sigs_up)
signaturesDown_df <- list_to_wide_df(flat_sigs_down)

xlsx_path <- file.path(save_path, "GR_signatures_up_down_31.10.2025.xlsx")

wb <- createWorkbook()
addWorksheet(wb, "signaturesUp")
addWorksheet(wb, "signaturesDown")
writeData(wb, "signaturesUp", signaturesUp_df)
writeData(wb, "signaturesDown", signaturesDown_df)
freezePane(wb, sheet = "signaturesUp", firstRow = TRUE)
freezePane(wb, sheet = "signaturesDown", firstRow = TRUE)
setColWidths(wb, sheet = "signaturesUp", cols = 1:ncol(signaturesUp_df), widths = "auto")
setColWidths(wb, sheet = "signaturesDown", cols = 1:ncol(signaturesDown_df), widths = "auto")
saveWorkbook(wb, xlsx_path, overwrite = TRUE)

cat("
✅ XLSX file created successfully:
   results_v2/GR_signatures/GR_signatures_31.10.2025/GR_signatures_up_down_31.10.2025.xlsx
   ├── Sheet 1: signaturesUp   (", length(flat_sigs_up),   " signatures)
   └── Sheet 2: signaturesDown (", length(flat_sigs_down), " signatures)
")
