# ============================================================
# 📘 Eksport list genów GR-zależnych (3 publikacje) do pliku XLSX
# ============================================================

library(openxlsx)

# --- 1️⃣ Przygotowanie list danych ---
sheet_names <- c(
  "NeuralSystemCells-up",
  "NeuralSystemCells-down",
  "BloodCells-up",
  "BloodCells-down",
  "PulmonaryCells-up",
  "PulmonaryCells-down"
)

# --- 2️⃣ Pobranie danych ---
gene_lists <- list(
  gene_summaryList_brain$n10_short_time$summary_up$enrichr_summary_source$genes_3pub %>% rev,
  gene_summaryList_brain$n10_short_time$summary_down$enrichr_summary_source$genes_3pub %>% rev,
  gene_summaryList_blood$n10_short_time$summary_up$enrichr_summary_source$genes_3pub %>% rev,
  gene_summaryList_blood$n10_short_time$summary_down$enrichr_summary_source$genes_3pub %>% rev,
  gene_summaryList_lung$n10_short_time$summary_up$enrichr_summary_source$genes_3pub %>% rev,
  gene_summaryList_lung$n10_short_time$summary_down$enrichr_summary_source$genes_3pub %>% rev
)

# --- 3️⃣ Utworzenie workbooka ---
wb <- createWorkbook()

# --- 4️⃣ Dodanie arkuszy i zapis danych ---
for (i in seq_along(gene_lists)) {
  df <- data.frame(hgnc_symbol = gene_lists[[i]])
  sheet_name <- sheet_names[i]
  
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet = sheet_name, x = df, headerStyle = createStyle(textDecoration = "bold"))
  
  # 🔒 Zablokuj pierwszy wiersz (nagłówek)
  freezePane(wb, sheet = sheet_name, firstRow = TRUE)
}

# --- 5️⃣ Zapis do pliku ---
saveWorkbook(wb, file = "results_v2/summary_geneLists_database/grSignatures/GRSignatures_genesPublication3_16.10.2025.xlsx", overwrite = TRUE)

message("✅ Zapisano plik: GR_signatures_genes3pub.xlsx")



# ============================================================
# 📘 Eksport list genów GR-zależnych (4 publikacje) do pliku XLSX
# ============================================================

library(openxlsx)

# --- 1️⃣ Przygotowanie list danych ---
sheet_names <- c(
  "NeuralSystemCells-up",
  "NeuralSystemCells-down",
  "BloodCells-up",
  "BloodCells-down",
  "PulmonaryCells-up",
  "PulmonaryCells-down"
)

# --- 2️⃣ Pobranie danych ---
gene_lists <- list(
  gene_summaryList_brain$n10_short_time$summary_up$enrichr_summary_source$genes_4pub %>% rev,
  gene_summaryList_brain$n10_short_time$summary_down$enrichr_summary_source$genes_4pub %>% rev,
  gene_summaryList_blood$n10_short_time$summary_up$enrichr_summary_source$genes_4pub %>% rev,
  gene_summaryList_blood$n10_short_time$summary_down$enrichr_summary_source$genes_4pub %>% rev,
  gene_summaryList_lung$n10_short_time$summary_up$enrichr_summary_source$genes_4pub %>% rev,
  gene_summaryList_lung$n10_short_time$summary_down$enrichr_summary_source$genes_4pub %>% rev
)

# --- 3️⃣ Utworzenie workbooka ---
wb <- createWorkbook()

# --- 4️⃣ Dodanie arkuszy i zapis danych ---
for (i in seq_along(gene_lists)) {
  df <- data.frame(hgnc_symbol = gene_lists[[i]])
  sheet_name <- sheet_names[i]
  
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet = sheet_name, x = df, headerStyle = createStyle(textDecoration = "bold"))
  
  # 🔒 Zablokuj pierwszy wiersz (nagłówek)
  freezePane(wb, sheet = sheet_name, firstRow = TRUE)
}

# --- 5️⃣ Zapis do pliku ---
saveWorkbook(wb, file = "results_v2/summary_geneLists_database/grSignatures/GRSignatures_genesPublication4_16.10.2025.xlsx", overwrite = TRUE)

message("✅ Zapisano plik: GRSignatures_genesPublication4_16.10.2025.xlsx")