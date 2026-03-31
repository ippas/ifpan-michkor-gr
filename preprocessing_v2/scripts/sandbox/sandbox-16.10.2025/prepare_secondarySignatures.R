gene_summaryList_brain$n10_short_time$summary_up$enrichr_summary_source$genes_3pub
gene_summaryList_blood$n10_short_time$summary_up$enrichr_summary_source$genes_3pub
gene_summaryList_lung$n10_short_time$summary_up$enrichr_summary_source$genes_3pub


list(
  gene_summaryList_brain$n10_short_time$summary_up$enrichr_summary_source$genes_3pub,
  gene_summaryList_blood$n10_short_time$summary_up$enrichr_summary_source$genes_3pub,
  gene_summaryList_lung$n10_short_time$summary_up$enrichr_summary_source$genes_3pub
) %>%
  unlist %>%
  table %>% as.data.frame() %>%
  filter(Freq > 1) -> repeated_up


list(
  gene_summaryList_brain$n10_short_time$summary_down$enrichr_summary_source$genes_3pub,
  gene_summaryList_blood$n10_short_time$summary_down$enrichr_summary_source$genes_3pub,
  gene_summaryList_lung$n10_short_time$summary_down$enrichr_summary_source$genes_3pub
) %>%
  unlist %>%
  table %>% as.data.frame() %>%
  filter(Freq > 1) -> repeated_down



repeatedTissueCells <- list(repeatedTissueCellsUp <- repeated_up,
                            repeatedTissueCellsDown <- repeated_down)



tissue_up_signatures_genes <- list(
  NeuralCellsUp = gene_summaryList_brain$n10_short_time$summary_up$enrichr_summary_source$genes_3pub,
  BloodCellsUp = gene_summaryList_blood$n10_short_time$summary_up$enrichr_summary_source$genes_3pub,
  LungCellsUp = gene_summaryList_lung$n10_short_time$summary_up$enrichr_summary_source$genes_3pub
)

# --- Wspólne geny (DOWN) ---
tissue_down_signatures_genes <- list(
  NeuralCellsDown = gene_summaryList_brain$n10_short_time$summary_down$enrichr_summary_source$genes_3pub,
  BloodCellsDown = gene_summaryList_blood$n10_short_time$summary_down$enrichr_summary_source$genes_3pub,
  LungCellsDown = gene_summaryList_lung$n10_short_time$summary_down$enrichr_summary_source$genes_3pub
)


# ============================================================
# 🔹 Usunięcie powtarzających się genów z sygnatur (UP)
# ============================================================

tissueCellsUpSignaturesMinusRepeated_genes <- list(
  NeuralCellsUp_minusRepeatedGenes = 
    tissue_up_signatures_genes$NeuralCellsUp[!tissue_up_signatures_genes$NeuralCellsUp %in% repeated_up],
  
  BloodCellsUp_minusRepeatedGenes = 
    tissue_up_signatures_genes$BloodCellsUp[!tissue_up_signatures_genes$BloodCellsUp %in% repeated_up],
  
  LungCellsUp_minusRepeatedGenes = 
    tissue_up_signatures_genes$LungCellsUp[!tissue_up_signatures_genes$LungCellsUp %in% repeated_up]
)



# ============================================================
# 🔹 Usunięcie powtarzających się genów z sygnatur (DOWN)
# ============================================================

tissueCellsDownSignaturesMinusRepeated_genes <- list(
  NeuralCellsDown_minusRepeatedGenes = 
    tissue_down_signatures_genes$NeuralCellsDown[!tissue_down_signatures_genes$NeuralCellsDown %in% repeated_down],
  
  BloodCellsDown_minusRepeatedGenes = 
    tissue_down_signatures_genes$BloodCellsDown[!tissue_down_signatures_genes$BloodCellsDown %in% repeated_down],
  
  LungCellsDown_minusRepeatedGenes = 
    tissue_down_signatures_genes$LungCellsDown[!tissue_down_signatures_genes$LungCellsDown %in% repeated_down]
)

signaturesMinusClusters_up <- list(
  
  minusClusterP = list(
    NeuralCellsUp_minusClusterP = 
      tissue_up_signatures_genes$NeuralCellsUp[!tissue_up_signatures_genes$NeuralCellsUp %in% marpiech_clusters$cluster_P],
    BloodCellsUp_minusClusterP = 
      tissue_up_signatures_genes$BloodCellsUp[!tissue_up_signatures_genes$BloodCellsUp %in% marpiech_clusters$cluster_P],
    LungCellsUp_minusClusterP = 
      tissue_up_signatures_genes$LungCellsUp[!tissue_up_signatures_genes$LungCellsUp %in% marpiech_clusters$cluster_P]
  ),
  
  minusClusterO = list(
    NeuralCellsUp_minusClusterO = 
      tissue_up_signatures_genes$NeuralCellsUp[!tissue_up_signatures_genes$NeuralCellsUp %in% marpiech_clusters$cluster_O],
    BloodCellsUp_minusClusterO = 
      tissue_up_signatures_genes$BloodCellsUp[!tissue_up_signatures_genes$BloodCellsUp %in% marpiech_clusters$cluster_O],
    LungCellsUp_minusClusterO = 
      tissue_up_signatures_genes$LungCellsUp[!tissue_up_signatures_genes$LungCellsUp %in% marpiech_clusters$cluster_O]
  ),
  
  minusClusterK = list(
    NeuralCellsUp_minusClusterK = 
      tissue_up_signatures_genes$NeuralCellsUp[!tissue_up_signatures_genes$NeuralCellsUp %in% marpiech_clusters$cluster_K],
    BloodCellsUp_minusClusterK = 
      tissue_up_signatures_genes$BloodCellsUp[!tissue_up_signatures_genes$BloodCellsUp %in% marpiech_clusters$cluster_K],
    LungCellsUp_minusClusterK = 
      tissue_up_signatures_genes$LungCellsUp[!tissue_up_signatures_genes$LungCellsUp %in% marpiech_clusters$cluster_K]
  ),
  
  minusClustersKPO = list(
    NeuralCellsUp_minusClustersKPO = 
      tissue_up_signatures_genes$NeuralCellsUp[!tissue_up_signatures_genes$NeuralCellsUp %in% unlist(marpiech_clusters[c("cluster_K", "cluster_P", "cluster_O")], use.names = FALSE)],
    BloodCellsUp_minusClustersKPO = 
      tissue_up_signatures_genes$BloodCellsUp[!tissue_up_signatures_genes$BloodCellsUp %in% unlist(marpiech_clusters[c("cluster_K", "cluster_P", "cluster_O")], use.names = FALSE)],
    LungCellsUp_minusClustersKPO = 
      tissue_up_signatures_genes$LungCellsUp[!tissue_up_signatures_genes$LungCellsUp %in% unlist(marpiech_clusters[c("cluster_K", "cluster_P", "cluster_O")], use.names = FALSE)]
  )
)



# ============================================================
# 🔹 Usunięcie genów należących do klastrów z sygnatur tkanek (DOWN)
# ============================================================

signaturesMinusClusters_down <- list(
  minusClusterD = list(
    NeuralCellsDown_minusClusterD = 
      tissue_down_signatures_genes$NeuralCellsDown[!tissue_down_signatures_genes$NeuralCellsDown %in% marpiech_clusters$cluster_D],
    BloodCellsDown_minusClusterD = 
      tissue_down_signatures_genes$BloodCellsDown[!tissue_down_signatures_genes$BloodCellsDown %in% marpiech_clusters$cluster_D],
    LungCellsDown_minusClusterD = 
      tissue_down_signatures_genes$LungCellsDown[!tissue_down_signatures_genes$LungCellsDown %in% marpiech_clusters$cluster_D]
  )
)

# ##############################################################################
# ---- prepare global signatures from full database ----
# ##############################################################################

# --- 2️⃣ Wstępne filtrowanie danych ---
filtered_data <- papers_data_preprocessing %>%
  filter(!(
    time %in% c(
      "240h",
      "720h",
      "240h_vs_720h",
      "9weeks",
      "3months",
      "10days",
      "14wekks",
      "8weeks",
      "168h",
      "240-528h",
      "14weeks",
      "3weeks",
      "7weeks",
      "672h",
      "main_effect_of_treatments"
    )
  )) %>%
  filter(!simple_tissue %in% c("Other", "placenta")) %>%
  filter(n_genes >= 10) %>%
  filter(regulation %in% c("up", "down")) %>% 
  mutate(simple_tissue_regulation = paste0(simple_tissue, "_", regulation))


# --- 🔹 Up-regulated genes ---
freq_table_up <- filtered_data %>%
  filter(regulation == "up") %>%
  group_by(simple_tissue_regulation) %>%
  nest() %>%
  mutate(genes = map(data, ~ unique(.x$hgnc_symbol))) %>%
  select(simple_tissue_regulation, genes) %>%
  pull(genes) %>%
  unlist() %>%
  table() %>%
  as.data.frame() %>%
  setNames(c("hgnc_symbol", "freq")) %>%
  count(freq, name = "n_genes") %>%
  arrange(desc(freq))

# --- 🔹 Down-regulated genes ---
freq_table_down <- filtered_data %>%
  filter(regulation == "down") %>%
  group_by(simple_tissue_regulation) %>%
  nest() %>%
  mutate(genes = map(data, ~ unique(.x$hgnc_symbol))) %>%
  select(simple_tissue_regulation, genes) %>%
  pull(genes) %>%
  unlist() %>%
  table() %>%
  as.data.frame() %>%
  setNames(c("hgnc_symbol", "freq")) %>%
  count(freq, name = "n_genes") %>%
  arrange(desc(freq))

# ======================================================
# 🥧 Pie chart: Up-regulated
# ======================================================
svg("results_v2/figures/summary_gr_database/tissuesGenes_piecharts/piechart_upArrangeTissuesGenes.svg",
    width = 5, height = 5, bg = "white")

ggplot(freq_table_up, aes(x = "", y = n_genes, fill = as.factor(freq))) +
  geom_bar(stat = "identity", color = "black", fill = "gray90") +
  coord_polar(theta = "y") +
  theme_void() +
  labs(title = "UP-regulated genes by frequency") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

dev.off()  # ✅ bardzo ważne: zamyka urządzenie graficzne

# ======================================================
# 🥧 Pie chart: Down-regulated
# ======================================================
svg("results_v2/figures/summary_gr_database/tissuesGenes_piecharts/piechart_downArrangeTissuesGenes.svg",
    width = 5, height = 5, bg = "white")

ggplot(freq_table_down, aes(x = "", y = n_genes, fill = as.factor(freq))) +
  geom_bar(stat = "identity", color = "black", fill = "gray90") +
  coord_polar(theta = "y") +
  theme_void() +
  labs(title = "DOWN-regulated genes by frequency") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

dev.off()  # zamyka plik SVG

# --- 3️⃣ Funkcja pomocnicza ---
get_genes_by_regulation_and_threshold <-
  function(df, regulation_type, threshold) {
    df %>%
      filter(regulation == regulation_type) %>%
      group_by(simple_tissue_regulation) %>%
      nest() %>%
      mutate(genes = map(data, ~ unique(.x$hgnc_symbol))) %>%
      select(-data) %>%
      pull(genes) %>%
      unlist() %>%
      table() %>%
      as.data.frame() %>%
      setNames(c("hgnc_symbol", "freq")) %>%
      filter(freq >= threshold) %>%
      arrange(desc(freq)) %>%
      pull(hgnc_symbol)
  }

# --- 4️⃣ Generowanie list genów globalnych ---
global_GR_genes <- list(
  globalUp4TissuesDerivedCells   = get_genes_by_regulation_and_threshold(filtered_data, "up",   4),
  globalUp5TissuesDerivedCells   = get_genes_by_regulation_and_threshold(filtered_data, "up",   5),
  globalUp6TissuesDerivedCells   = get_genes_by_regulation_and_threshold(filtered_data, "up",   6),
  globalDown4TissuesDerivedCells = get_genes_by_regulation_and_threshold(filtered_data, "down", 4),
  globalDown5TissuesDerivedCells = get_genes_by_regulation_and_threshold(filtered_data, "down", 5),
  globalDown6TissuesDerivedCells = get_genes_by_regulation_and_threshold(filtered_data, "down", 6)
)


# --- 5️⃣ Podsumowanie liczby genów ---
cat("\n📊 Liczba genów w każdej kategorii:\n")
print(sapply(general_GR_genes, length))


signaturesMinusGlobal <- list(
  
  # --- Globalne dla 4 tkanek ---
  minusGlobalUp4TissuesDerivedCells = list(
    NeuralCellsUp_minusGlobalUp4TissuesDerivedCells = 
      tissue_up_signatures_genes$NeuralCellsUp[!tissue_up_signatures_genes$NeuralCellsUp %in% global_GR_genes$globalUp_4tissuesDerivedCells],
    BloodCellsUp_minusGlobalUp4TissuesDerivedCells = 
      tissue_up_signatures_genes$BloodCellsUp[!tissue_up_signatures_genes$BloodCellsUp %in% global_GR_genes$globalUp_4tissuesDerivedCells],
    LungCellsUp_minusGlobalUp4TissuesDerivedCells = 
      tissue_up_signatures_genes$LungCellsUp[!tissue_up_signatures_genes$LungCellsUp %in% global_GR_genes$globalUp_4tissuesDerivedCells]
  ),
  
  minusGlobalDown4TissuesDerivedCells = list(
    NeuralCellsDown_minusGlobalDown4TissuesDerivedCells = 
      tissue_down_signatures_genes$NeuralCellsDown[!tissue_down_signatures_genes$NeuralCellsDown %in% global_GR_genes$globalDown_4tissuesDerivedCells],
    BloodCellsDown_minusGlobalDown4TissuesDerivedCells = 
      tissue_down_signatures_genes$BloodCellsDown[!tissue_down_signatures_genes$BloodCellsDown %in% global_GR_genes$globalDown_4tissuesDerivedCells],
    LungCellsDown_minusGlobalDown4TissuesDerivedCells = 
      tissue_down_signatures_genes$LungCellsDown[!tissue_down_signatures_genes$LungCellsDown %in% global_GR_genes$globalDown_4tissuesDerivedCells]
  ),
  
  # --- Globalne dla 5 tkanek ---
  minusGlobalUp5TissuesDerivedCells = list(
    NeuralCellsUp_minusGlobalUp5TissuesDerivedCells = 
      tissue_up_signatures_genes$NeuralCellsUp[!tissue_up_signatures_genes$NeuralCellsUp %in% global_GR_genes$globalUp_5tissuesDerivedCells],
    BloodCellsUp_minusGlobalUp5TissuesDerivedCells = 
      tissue_up_signatures_genes$BloodCellsUp[!tissue_up_signatures_genes$BloodCellsUp %in% global_GR_genes$globalUp_5tissuesDerivedCells],
    LungCellsUp_minusGlobalUp5TissuesDerivedCells = 
      tissue_up_signatures_genes$LungCellsUp[!tissue_up_signatures_genes$LungCellsUp %in% global_GR_genes$globalUp_5tissuesDerivedCells]
  ),
  
  minusGlobalDown5TissuesDerivedCells = list(
    NeuralCellsDown_minusGlobalDown5TissuesDerivedCells = 
      tissue_down_signatures_genes$NeuralCellsDown[!tissue_down_signatures_genes$NeuralCellsDown %in% global_GR_genes$globalDown_5tissuesDerivedCells],
    BloodCellsDown_minusGlobalDown5TissuesDerivedCells = 
      tissue_down_signatures_genes$BloodCellsDown[!tissue_down_signatures_genes$BloodCellsDown %in% global_GR_genes$globalDown_5tissuesDerivedCells],
    LungCellsDown_minusGlobalDown5TissuesDerivedCells = 
      tissue_down_signatures_genes$LungCellsDown[!tissue_down_signatures_genes$LungCellsDown %in% global_GR_genes$globalDown_5tissuesDerivedCells]
  ),
  
  # --- Globalne dla 6 tkanek ---
  minusGlobalUp6TissuesDerivedCells = list(
    NeuralCellsUp_minusGlobalUp6TissuesDerivedCells = 
      tissue_up_signatures_genes$NeuralCellsUp[!tissue_up_signatures_genes$NeuralCellsUp %in% global_GR_genes$globalUp_6tissuesDerivedCells],
    BloodCellsUp_minusGlobalUp6TissuesDerivedCells = 
      tissue_up_signatures_genes$BloodCellsUp[!tissue_up_signatures_genes$BloodCellsUp %in% global_GR_genes$globalUp_6tissuesDerivedCells],
    LungCellsUp_minusGlobalUp6TissuesDerivedCells = 
      tissue_up_signatures_genes$LungCellsUp[!tissue_up_signatures_genes$LungCellsUp %in% global_GR_genes$globalUp_6tissuesDerivedCells]
  ),
  
  minusGlobalDown6TissuesDerivedCells = list(
    NeuralCellsDown_minusGlobalDown6TissuesDerivedCells = 
      tissue_down_signatures_genes$NeuralCellsDown[!tissue_down_signatures_genes$NeuralCellsDown %in% global_GR_genes$globalDown_6tissuesDerivedCells],
    BloodCellsDown_minusGlobalDown6TissuesDerivedCells = 
      tissue_down_signatures_genes$BloodCellsDown[!tissue_down_signatures_genes$BloodCellsDown %in% global_GR_genes$globalDown_6tissuesDerivedCells],
    LungCellsDown_minusGlobalDown6TissuesDerivedCells = 
      tissue_down_signatures_genes$LungCellsDown[!tissue_down_signatures_genes$LungCellsDown %in% global_GR_genes$globalDown_6tissuesDerivedCells]
  )
)

allGrSignatures_17.10.2025 <- list(
  tissue_up_signatures_genes     = tissue_up_signatures_genes,
  tissue_down_signatures_genes   = tissue_down_signatures_genes,
  signaturesMinusClusters_up     = signaturesMinusClusters_up,
  signaturesMinusClusters_down   = signaturesMinusClusters_down,
  global_GR_genes                = global_GR_genes,
  signaturesMinusGlobal          = signaturesMinusGlobal,
  tissueCellsUpSignaturesMinusRepeated_genes = tissueCellsUpSignaturesMinusRepeated_genes,
  tissueCellsDownSignaturesMinusRepeated_genes = tissueCellsDownSignaturesMinusRepeated_genes,
  repeatedTissueCells = repeatedTissueCells 
)

flatten_named_list <- function(x, parent_name = NULL) {
  out <- list()
  
  if (is.list(x)) {
    for (nm in names(x)) {
      # Tworzymy pełną nazwę ścieżki
      full_name <- if (is.null(parent_name)) nm else paste(parent_name, nm, sep = ".")
      # Rekurencja
      out <- c(out, flatten_named_list(x[[nm]], full_name))
    }
  } else if (is.atomic(x)) {
    # Ostatni człon nazwy po kropce
    short_name <- if (!is.null(parent_name)) sub(".*\\.", "", parent_name) else parent_name
    out[[short_name]] <- as.character(x)
  }
  
  return(out)
}
# ============================================================
# 🧬 Zastosowanie do allGrSignatures_17.10.2025
# ============================================================

flat_allGrSignatures_17.10.2025 <- flatten_named_list(allGrSignatures_17.10.2025)


df_allGrSignatures_17.10.2025 <- map2_df(
  flat_allGrSignatures_17.10.2025,
  names(flat_allGrSignatures_17.10.2025),
  ~ tibble(
    hgnc_symbol = .x,
    signature_name = .y
  )
)


# --- 6️⃣ Porządkowanie środowiska: zostaw tylko wynikową listę ---
rm(filtered_data, get_genes_by_regulation_and_threshold)
