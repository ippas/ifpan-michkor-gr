# ============================================================
# 🧬 Preparation and cleanup of GR-dependent gene signatures
# ============================================================

library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

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
# 7️⃣ Combine all objects into one structured list
# ============================================================

allGrSignatures_17.10.2025 <- list(
  tissue_up_signatures_genes = tissue_up_signatures_genes,
  tissue_down_signatures_genes = tissue_down_signatures_genes,
  tissueCellsUpSignaturesMinusRepeated_genes = tissueCellsUpSignaturesMinusRepeated_genes,
  tissueCellsDownSignaturesMinusRepeated_genes = tissueCellsDownSignaturesMinusRepeated_genes,
  repeatedTissueCells = repeatedTissueCells,
  signaturesMinusClusters_up = signaturesMinusClusters_up,
  signaturesMinusClusters_down = signaturesMinusClusters_down,
  global_GR_genes = global_GR_genes,
  signaturesMinusGlobal = signaturesMinusGlobal
)

# ============================================================
# 8️⃣ Flatten nested structure into a data frame
# ============================================================

flatten_named_list <- function(x, parent = character(), sep = ".") {
  if (is.list(x)) {
    # rekurencja po elementach listy z dołączaniem nazw
    purrr::imap(x, ~ flatten_named_list(.x, c(parent, .y), sep)) |> unlist(recursive = FALSE)
  } else if (is.atomic(x)) {
    nm <- paste(parent, collapse = sep)            # pełna ścieżka, np. "signaturesMinusGlobal.minusGlobalUp4TissuesDerivedCells.NeuralCellsUp_minusGlobalUp4TissuesDerivedCells"
    setNames(list(as.character(x)), nm)
  } else {
    list()
  }
}

# flat_allGrSignatures_17.10.2025 <- flatten_named_list(allGrSignatures_17.10.2025)

names(flat_allGrSignatures_17.10.2025) <- sapply(names(flat_allGrSignatures_17.10.2025), function(nm) {
  parts <- unlist(strsplit(nm, "\\."))
  new_name <- if (length(parts) >= 2) {
    paste0(parts[length(parts) - 1], "_", parts[length(parts)])
  } else {
    nm
  }
  
  # remove unwanted prefixes
  new_name <- gsub("^tissue_up_signatures_genes_", "", new_name)
  new_name <- gsub("^tissue_down_signatures_genes_", "", new_name)
  
  return(new_name)
})

df_allGrSignatures_17.10.2025 <- map2_df(
  flat_allGrSignatures_17.10.2025,
  names(flat_allGrSignatures_17.10.2025),
  \(x, nm) tibble(hgnc_symbol = x, signature_name = nm)
)

# ============================================================
# ✅ Final environment cleanup (keep only key results)
# ============================================================

rm(
  repeated_up, repeated_down, filtered_data,
  remove_repeated, remove_cluster_genes,
  remove_global_genes, get_genes_by_regulation_and_threshold
)



# ============================================================
# 💾 Save GR signature objects (RDS, TSV, XLSX)
# ============================================================

library(openxlsx)
library(readr)

# --- 1️⃣ Define save path ---
save_path <- "results_v2/GR_signatures/GR_signatures_17.10.2025"

# --- 2️⃣ Save RDS objects (same names as in environment) ---
saveRDS(flat_allGrSignatures_17.10.2025,
        file = file.path(save_path, "flat_allGrSignatures_17.10.2025.rds"))

saveRDS(allGrSignatures_17.10.2025,
        file = file.path(save_path, "allGrSignatures_17.10.2025.rds"))

saveRDS(df_allGrSignatures_17.10.2025,
        file = file.path(save_path, "df_allGrSignatures_17.10.2025.rds"))

# --- 3️⃣ Export df_allGrSignatures_17.10.2025 as TSV ---
write_tsv(df_allGrSignatures_17.10.2025,
          file.path(save_path, "df_allGrSignatures_17.10.2025.tsv"))

# --- 4️⃣ Export df_allGrSignatures_17.10.2025 as XLSX ---
xlsx_file <- file.path(save_path, "df_allGrSignatures_17.10.2025.xlsx")
wb <- createWorkbook()
addWorksheet(wb, "GR_signatures")
writeData(wb, "GR_signatures", df_allGrSignatures_17.10.2025)
saveWorkbook(wb, xlsx_file, overwrite = TRUE)

# --- 5️⃣ Confirmation message ---
cat("
✅ Saved GR signature objects to: results_v2/GR_signatures/GR_signatures_17.10.2025/
   ├── flat_allGrSignatures_17.10.2025.rds
   ├── allGrSignatures_17.10.2025.rds
   ├── df_allGrSignatures_17.10.2025.rds
   ├── df_allGrSignatures_17.10.2025.tsv
   └── df_allGrSignatures_17.10.2025.xlsx
")



# ============================================================
# 📘 Export UP/DOWN signatures from flat_allGrSignatures_17.10.2025 to XLSX
# ============================================================

library(openxlsx)
library(purrr)

# --- 1️⃣ Select UP and DOWN signatures directly ---
flat_sigs_up   <- flat_allGrSignatures_17.10.2025[grepl("Up",   names(flat_allGrSignatures_17.10.2025))]
flat_sigs_down <- flat_allGrSignatures_17.10.2025[grepl("Down", names(flat_allGrSignatures_17.10.2025))]

# --- 2️⃣ Helper: convert list of vectors to rectangular data frame ---
list_to_wide_df <- function(lst) {
  if (length(lst) == 0) return(data.frame())
  max_len <- max(lengths(lst))
  data.frame(lapply(lst, function(x) c(x, rep(NA, max_len - length(x)))), check.names = FALSE)
}

signaturesUp_df   <- list_to_wide_df(flat_sigs_up)
signaturesDown_df <- list_to_wide_df(flat_sigs_down)

# --- 3️⃣ Save to XLSX ---
xlsx_path <- "results_v2/GR_signatures/GR_signatures_17.10.2025/GR_signatures_up_down_17.10.2025.xlsx"

wb <- createWorkbook()
addWorksheet(wb, "signaturesUp")
addWorksheet(wb, "signaturesDown")

writeData(wb, "signaturesUp",   signaturesUp_df)
writeData(wb, "signaturesDown", signaturesDown_df)

# --- 🔒 Freeze the first row in both sheets ---
freezePane(wb, sheet = "signaturesUp", firstRow = TRUE)
freezePane(wb, sheet = "signaturesDown", firstRow = TRUE)

# --- 📏 Auto-adjust column widths to longest entry ---
setColWidths(wb, sheet = "signaturesUp",   cols = 1:ncol(signaturesUp_df),   widths = "auto")
setColWidths(wb, sheet = "signaturesDown", cols = 1:ncol(signaturesDown_df), widths = "auto")

saveWorkbook(wb, xlsx_path, overwrite = TRUE)

# --- 4️⃣ Confirmation ---
cat("
✅ XLSX file created successfully:
   results_v2/GR_signatures/GR_signatures_17.10.2025/GR_signatures_up_down_17.10.2025.xlsx
   ├── Sheet 1: signaturesUp   (", length(flat_sigs_up),   " signatures)
   └── Sheet 2: signaturesDown (", length(flat_sigs_down), " signatures)
")


rm(
  flat_sigs_up,
  flat_sigs_down,
  list_to_wide_df,
  signaturesUp_df,
  signaturesDown_df,
  wb,
  xlsx_path
)
