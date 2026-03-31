# 🧠 „Neural system cells”, 🩸 „Blood cells”, 🫁 „Pulmonary cells”

# ##############################################################################
# ---- CellMarker_2024 ----
# ##############################################################################

library(dplyr)
library(openxlsx)

# --- 1️⃣ Przygotowanie filtrów ---
filter_results <- function(df) {
  df %>%
    filter(n_genes >= 3) %>%
    filter(Adjusted.P.value < 0.05)
}

# --- 2️⃣ Zestawy danych ---
NeuralSystemCells_up   <- filter_results(gene_summaryList_brain$n10_short_time$summary_up$enrichr_summary_source$enrichr_3pub$CellMarker_2024)
NeuralSystemCells_down <- filter_results(gene_summaryList_brain$n10_short_time$summary_down$enrichr_summary_source$enrichr_3pub$CellMarker_2024)

BloodCells_up   <- filter_results(gene_summaryList_blood$n10_short_time$summary_up$enrichr_summary_source$enrichr_3pub$CellMarker_2024)
BloodCells_down <- filter_results(gene_summaryList_blood$n10_short_time$summary_down$enrichr_summary_source$enrichr_3pub$CellMarker_2024)

PulmonaryCells_up   <- filter_results(gene_summaryList_lung$n10_short_time$summary_up$enrichr_summary_source$enrichr_3pub$CellMarker_2024)
PulmonaryCells_down <- filter_results(gene_summaryList_lung$n10_short_time$summary_down$enrichr_summary_source$enrichr_3pub$CellMarker_2024)

# --- 3️⃣ Tworzenie workbooka ---
wb <- createWorkbook()

sheet_names <- c(
  "NeuralSystemCells-up",
  "NeuralSystemCells-down",
  "BloodCells-up",
  "BloodCells-down",
  "PulmonaryCells-up",
  "PulmonaryCells-down"
)

datasets <- list(
  NeuralSystemCells_up,
  NeuralSystemCells_down,
  BloodCells_up,
  BloodCells_down,
  PulmonaryCells_up,
  PulmonaryCells_down
)

# --- 4️⃣ Dodanie arkuszy i wpisanie danych ---
for (i in seq_along(sheet_names)) {
  addWorksheet(wb, sheet_names[i])
  writeData(wb, sheet_names[i], datasets[[i]])
  freezePane(wb, sheet = sheet_names[i], firstRow = TRUE)  # 🔒 zablokowany nagłówek
  setColWidths(wb, sheet = sheet_names[i], cols = 1:ncol(datasets[[i]]), widths = "auto")  # 📏 auto szerokość
}

# --- 5️⃣ Ścieżka i zapis ---
output_dir <- "results_v2/summary_geneLists_database/functional_validation_grSignatures"
output_file <- file.path(output_dir, "CellMarker2024_functionalValidationSignaturesPublication3_fdr0.05overlap3_16.10.2025.xlsx")

# Utwórz folder jeśli nie istnieje
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

saveWorkbook(wb, file = output_file, overwrite = TRUE)

cat("✅ Zapisano do pliku:\n", output_file, "\n(nagłówki zablokowane, kolumny dopasowane)\n")

# ##############################################################################
# ---- ChEA 2022 ----
# ##############################################################################

# --- 2️⃣ Zestawy danych ---
NeuralSystemCells_up   <- filter_results(gene_summaryList_brain$n10_short_time$summary_up$enrichr_summary_source$enrichr_3pub$ChEA_2022)
NeuralSystemCells_down <- filter_results(gene_summaryList_brain$n10_short_time$summary_down$enrichr_summary_source$enrichr_3pub$ChEA_2022)

BloodCells_up   <- filter_results(gene_summaryList_blood$n10_short_time$summary_up$enrichr_summary_source$enrichr_3pub$ChEA_2022)
BloodCells_down <- filter_results(gene_summaryList_blood$n10_short_time$summary_down$enrichr_summary_source$enrichr_3pub$ChEA_2022)

PulmonaryCells_up   <- filter_results(gene_summaryList_lung$n10_short_time$summary_up$enrichr_summary_source$enrichr_3pub$ChEA_2022)
PulmonaryCells_down <- filter_results(gene_summaryList_lung$n10_short_time$summary_down$enrichr_summary_source$enrichr_3pub$ChEA_2022)

# --- 3️⃣ Tworzenie workbooka ---
wb <- createWorkbook()

sheet_names <- c(
  "NeuralSystemCells-up",
  "NeuralSystemCells-down",
  "BloodCells-up",
  "BloodCells-down",
  "PulmonaryCells-up",
  "PulmonaryCells-down"
)

datasets <- list(
  NeuralSystemCells_up,
  NeuralSystemCells_down,
  BloodCells_up,
  BloodCells_down,
  PulmonaryCells_up,
  PulmonaryCells_down
)

# --- 4️⃣ Dodanie arkuszy i wpisanie danych ---
for (i in seq_along(sheet_names)) {
  addWorksheet(wb, sheet_names[i])
  writeData(wb, sheet_names[i], datasets[[i]])
  freezePane(wb, sheet = sheet_names[i], firstRow = TRUE)  # 🔒 zablokowany nagłówek
  setColWidths(wb, sheet = sheet_names[i], cols = 1:ncol(datasets[[i]]), widths = "auto")  # 📏 auto szerokość
}

# --- 5️⃣ Ścieżka i zapis ---
output_dir <- "results_v2/summary_geneLists_database/functional_validation_grSignatures"
output_file <- file.path(output_dir, "ChEA2022_functionalValidationSignaturesPublication3_fdr0.05overlap3_16.10.2025.xlsx")

# Utwórz folder jeśli nie istnieje
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

saveWorkbook(wb, file = output_file, overwrite = TRUE)

cat("✅ Zapisano do pliku:\n", output_file, "\n(nagłówki zablokowane, kolumny dopasowane)\n")


# ##############################################################################
# ---- LINCS_L1000_Chem_Pert_Consensus_Sigs2 ----
# ##############################################################################

# --- 2️⃣ Zestawy danych ---
NeuralSystemCells_up   <- filter_results(gene_summaryList_brain$n10_short_time$summary_up$enrichr_summary_source$enrichr_3pub$LINCS_L1000_Chem_Pert_Consensus_Sigs)
NeuralSystemCells_down <- filter_results(gene_summaryList_brain$n10_short_time$summary_down$enrichr_summary_source$enrichr_3pub$LINCS_L1000_Chem_Pert_Consensus_Sigs)

BloodCells_up   <- filter_results(gene_summaryList_blood$n10_short_time$summary_up$enrichr_summary_source$enrichr_3pub$LINCS_L1000_Chem_Pert_Consensus_Sigs)
BloodCells_down <- filter_results(gene_summaryList_blood$n10_short_time$summary_down$enrichr_summary_source$enrichr_3pub$LINCS_L1000_Chem_Pert_Consensus_Sigs)

PulmonaryCells_up   <- filter_results(gene_summaryList_lung$n10_short_time$summary_up$enrichr_summary_source$enrichr_3pub$LINCS_L1000_Chem_Pert_Consensus_Sigs)
PulmonaryCells_down <- filter_results(gene_summaryList_lung$n10_short_time$summary_down$enrichr_summary_source$enrichr_3pub$LINCS_L1000_Chem_Pert_Consensus_Sigs)

# --- 3️⃣ Tworzenie workbooka ---
wb <- createWorkbook()

sheet_names <- c(
  "NeuralSystemCells-up",
  "NeuralSystemCells-down",
  "BloodCells-up",
  "BloodCells-down",
  "PulmonaryCells-up",
  "PulmonaryCells-down"
)

datasets <- list(
  NeuralSystemCells_up,
  NeuralSystemCells_down,
  BloodCells_up,
  BloodCells_down,
  PulmonaryCells_up,
  PulmonaryCells_down
)

# --- 4️⃣ Dodanie arkuszy i wpisanie danych ---
for (i in seq_along(sheet_names)) {
  addWorksheet(wb, sheet_names[i])
  writeData(wb, sheet_names[i], datasets[[i]])
  freezePane(wb, sheet = sheet_names[i], firstRow = TRUE)  # 🔒 zablokowany nagłówek
  setColWidths(wb, sheet = sheet_names[i], cols = 1:ncol(datasets[[i]]), widths = "auto")  # 📏 auto szerokość
}

# --- 5️⃣ Ścieżka i zapis ---
output_dir <- "results_v2/summary_geneLists_database/functional_validation_grSignatures"
output_file <- file.path(output_dir, "LINCSL1000_functionalValidationSignaturesPublication3_fdr0.05overlap3_16.10.2025.xlsx")

# Utwórz folder jeśli nie istnieje
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

saveWorkbook(wb, file = output_file, overwrite = TRUE)

cat("✅ Zapisano do pliku:\n", output_file, "\n(nagłówki zablokowane, kolumny dopasowane)\n")

# clean environments
rm(
  NeuralSystemCells_up,
  NeuralSystemCells_down,
  BloodCells_up,
  BloodCells_down,
  PulmonaryCells_up,
  PulmonaryCells_down,
  sheet_names,
  datasets,
  wb,
  output_dir,
  output_file,
  filter_results
)
gc()  # 🧹 zwolnienie pamięci


# ##############################################################################
# ---- for 4 publication ----
# ##############################################################################

NeuralSystemCells_up   <- filter_results(gene_summaryList_brain$n10_short_time$summary_up$enrichr_summary_source$enrichr_4pub$CellMarker_2024)
NeuralSystemCells_down <- filter_results(gene_summaryList_brain$n10_short_time$summary_down$enrichr_summary_source$enrichr_4pub$CellMarker_2024)

BloodCells_up   <- filter_results(gene_summaryList_blood$n10_short_time$summary_up$enrichr_summary_source$enrichr_4pub$CellMarker_2024)
BloodCells_down <- filter_results(gene_summaryList_blood$n10_short_time$summary_down$enrichr_summary_source$enrichr_4pub$CellMarker_2024)

PulmonaryCells_up   <- filter_results(gene_summaryList_lung$n10_short_time$summary_up$enrichr_summary_source$enrichr_4pub$CellMarker_2024)
PulmonaryCells_down <- filter_results(gene_summaryList_lung$n10_short_time$summary_down$enrichr_summary_source$enrichr_4pub$CellMarker_2024)

wb <- createWorkbook()
sheet_names <- c(
  "NeuralSystemCells-up", "NeuralSystemCells-down",
  "BloodCells-up", "BloodCells-down",
  "PulmonaryCells-up", "PulmonaryCells-down"
)
datasets <- list(
  NeuralSystemCells_up, NeuralSystemCells_down,
  BloodCells_up, BloodCells_down,
  PulmonaryCells_up, PulmonaryCells_down
)

for (i in seq_along(sheet_names)) {
  addWorksheet(wb, sheet_names[i])
  writeData(wb, sheet_names[i], datasets[[i]])
  freezePane(wb, sheet = sheet_names[i], firstRow = TRUE)
  setColWidths(wb, sheet = sheet_names[i], cols = 1:ncol(datasets[[i]]), widths = "auto")
}

output_dir <- "results_v2/summary_geneLists_database/functional_validation_grSignatures"
output_file <- file.path(output_dir, "CellMarker2024_functionalValidationSignaturesPublication4_fdr0.05overlap3_16.10.2025.xlsx")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
saveWorkbook(wb, file = output_file, overwrite = TRUE)
cat("✅ Zapisano:", output_file, "\n")

# ============================================================
# ---- ChEA_2022 ----
# ============================================================

NeuralSystemCells_up   <- filter_results(gene_summaryList_brain$n10_short_time$summary_up$enrichr_summary_source$enrichr_4pub$ChEA_2022)
NeuralSystemCells_down <- filter_results(gene_summaryList_brain$n10_short_time$summary_down$enrichr_summary_source$enrichr_4pub$ChEA_2022)

BloodCells_up   <- filter_results(gene_summaryList_blood$n10_short_time$summary_up$enrichr_summary_source$enrichr_4pub$ChEA_2022)
BloodCells_down <- filter_results(gene_summaryList_blood$n10_short_time$summary_down$enrichr_summary_source$enrichr_4pub$ChEA_2022)

PulmonaryCells_up   <- filter_results(gene_summaryList_lung$n10_short_time$summary_up$enrichr_summary_source$enrichr_4pub$ChEA_2022)
PulmonaryCells_down <- filter_results(gene_summaryList_lung$n10_short_time$summary_down$enrichr_summary_source$enrichr_4pub$ChEA_2022)

wb <- createWorkbook()
datasets <- list(
  NeuralSystemCells_up, NeuralSystemCells_down,
  BloodCells_up, BloodCells_down,
  PulmonaryCells_up, PulmonaryCells_down
)

for (i in seq_along(sheet_names)) {
  addWorksheet(wb, sheet_names[i])
  writeData(wb, sheet_names[i], datasets[[i]])
  freezePane(wb, sheet = sheet_names[i], firstRow = TRUE)
  setColWidths(wb, sheet = sheet_names[i], cols = 1:ncol(datasets[[i]]), widths = "auto")
}

output_file <- file.path(output_dir, "ChEA2022_functionalValidationSignaturesPublication4_fdr0.05overlap3_16.10.2025.xlsx")
saveWorkbook(wb, file = output_file, overwrite = TRUE)
cat("✅ Zapisano:", output_file, "\n")

# ============================================================
# ---- LINCS_L1000_Chem_Pert_Consensus_Sigs ----
# ============================================================

NeuralSystemCells_up   <- filter_results(gene_summaryList_brain$n10_short_time$summary_up$enrichr_summary_source$enrichr_4pub$LINCS_L1000_Chem_Pert_Consensus_Sigs)
NeuralSystemCells_down <- filter_results(gene_summaryList_brain$n10_short_time$summary_down$enrichr_summary_source$enrichr_4pub$LINCS_L1000_Chem_Pert_Consensus_Sigs)

BloodCells_up   <- filter_results(gene_summaryList_blood$n10_short_time$summary_up$enrichr_summary_source$enrichr_4pub$LINCS_L1000_Chem_Pert_Consensus_Sigs)
BloodCells_down <- filter_results(gene_summaryList_blood$n10_short_time$summary_down$enrichr_summary_source$enrichr_4pub$LINCS_L1000_Chem_Pert_Consensus_Sigs)

PulmonaryCells_up   <- filter_results(gene_summaryList_lung$n10_short_time$summary_up$enrichr_summary_source$enrichr_4pub$LINCS_L1000_Chem_Pert_Consensus_Sigs)
PulmonaryCells_down <- filter_results(gene_summaryList_lung$n10_short_time$summary_down$enrichr_summary_source$enrichr_4pub$LINCS_L1000_Chem_Pert_Consensus_Sigs)

wb <- createWorkbook()
datasets <- list(
  NeuralSystemCells_up, NeuralSystemCells_down,
  BloodCells_up, BloodCells_down,
  PulmonaryCells_up, PulmonaryCells_down
)

for (i in seq_along(sheet_names)) {
  addWorksheet(wb, sheet_names[i])
  writeData(wb, sheet_names[i], datasets[[i]])
  freezePane(wb, sheet = sheet_names[i], firstRow = TRUE)
  setColWidths(wb, sheet = sheet_names[i], cols = 1:ncol(datasets[[i]]), widths = "auto")
}

output_file <- file.path(output_dir, "LINCSL1000_functionalValidationSignaturesPublication4_fdr0.05overlap3_16.10.2025.xlsx")
saveWorkbook(wb, file = output_file, overwrite = TRUE)
cat("✅ Zapisano:", output_file, "\n")

# ============================================================
# 🧹 Czyszczenie środowiska
# ============================================================

rm(
  NeuralSystemCells_up, NeuralSystemCells_down,
  BloodCells_up, BloodCells_down,
  PulmonaryCells_up, PulmonaryCells_down,
  sheet_names, datasets, wb,
  output_dir, output_file, filter_results
)
gc()  # zwolnienie pamięci

cat("🧹 Środowisko wyczyszczone, wszystko zapisane poprawnie.\n")

# ##############################################################################
# ---- for 5 publication ----
# ##############################################################################
# --- 1️⃣ Przygotowanie filtrów ---
filter_results <- function(df) {
  df %>%
    filter(n_genes >= 3) %>%
    filter(Adjusted.P.value < 0.05)
}

# ============================================================
# ---- CellMarker_2024 ----
# ============================================================

NeuralSystemCells_up   <- filter_results(gene_summaryList_brain$n10_short_time$summary_up$enrichr_summary_source$enrichr_5pub$CellMarker_2024)
NeuralSystemCells_down <- filter_results(gene_summaryList_brain$n10_short_time$summary_down$enrichr_summary_source$enrichr_5pub$CellMarker_2024)

BloodCells_up   <- filter_results(gene_summaryList_blood$n10_short_time$summary_up$enrichr_summary_source$enrichr_5pub$CellMarker_2024)
BloodCells_down <- filter_results(gene_summaryList_blood$n10_short_time$summary_down$enrichr_summary_source$enrichr_5pub$CellMarker_2024)

PulmonaryCells_up   <- filter_results(gene_summaryList_lung$n10_short_time$summary_up$enrichr_summary_source$enrichr_5pub$CellMarker_2024)
PulmonaryCells_down <- filter_results(gene_summaryList_lung$n10_short_time$summary_down$enrichr_summary_source$enrichr_5pub$CellMarker_2024)

wb <- createWorkbook()
sheet_names <- c(
  "NeuralSystemCells-up", "NeuralSystemCells-down",
  "BloodCells-up", "BloodCells-down",
  "PulmonaryCells-up", "PulmonaryCells-down"
)
datasets <- list(
  NeuralSystemCells_up, NeuralSystemCells_down,
  BloodCells_up, BloodCells_down,
  PulmonaryCells_up, PulmonaryCells_down
)

for (i in seq_along(sheet_names)) {
  addWorksheet(wb, sheet_names[i])
  writeData(wb, sheet_names[i], datasets[[i]])
  freezePane(wb, sheet = sheet_names[i], firstRow = TRUE)
  setColWidths(wb, sheet = sheet_names[i], cols = 1:ncol(datasets[[i]]), widths = "auto")
}

output_dir <- "results_v2/summary_geneLists_database/functional_validation_grSignatures"
output_file <- file.path(output_dir, "CellMarker2024_functionalValidationSignaturesPublication5_fdr0.05overlap3_16.10.2025.xlsx")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
saveWorkbook(wb, file = output_file, overwrite = TRUE)
cat("✅ Zapisano:", output_file, "\n")

# ============================================================
# ---- ChEA_2022 ----
# ============================================================

NeuralSystemCells_up   <- filter_results(gene_summaryList_brain$n10_short_time$summary_up$enrichr_summary_source$enrichr_5pub$ChEA_2022)
NeuralSystemCells_down <- filter_results(gene_summaryList_brain$n10_short_time$summary_down$enrichr_summary_source$enrichr_5pub$ChEA_2022)

BloodCells_up   <- filter_results(gene_summaryList_blood$n10_short_time$summary_up$enrichr_summary_source$enrichr_5pub$ChEA_2022)
BloodCells_down <- filter_results(gene_summaryList_blood$n10_short_time$summary_down$enrichr_summary_source$enrichr_5pub$ChEA_2022)

PulmonaryCells_up   <- filter_results(gene_summaryList_lung$n10_short_time$summary_up$enrichr_summary_source$enrichr_5pub$ChEA_2022)
PulmonaryCells_down <- filter_results(gene_summaryList_lung$n10_short_time$summary_down$enrichr_summary_source$enrichr_5pub$ChEA_2022)

wb <- createWorkbook()
datasets <- list(
  NeuralSystemCells_up, NeuralSystemCells_down,
  BloodCells_up, BloodCells_down,
  PulmonaryCells_up, PulmonaryCells_down
)

for (i in seq_along(sheet_names)) {
  addWorksheet(wb, sheet_names[i])
  writeData(wb, sheet_names[i], datasets[[i]])
  freezePane(wb, sheet = sheet_names[i], firstRow = TRUE)
  setColWidths(wb, sheet = sheet_names[i], cols = 1:ncol(datasets[[i]]), widths = "auto")
}

output_file <- file.path(output_dir, "ChEA2022_functionalValidationSignaturesPublication5_fdr0.05overlap3_16.10.2025.xlsx")
saveWorkbook(wb, file = output_file, overwrite = TRUE)
cat("✅ Zapisano:", output_file, "\n")

# ============================================================
# ---- LINCS_L1000_Chem_Pert_Consensus_Sigs ----
# ============================================================

NeuralSystemCells_up   <- filter_results(gene_summaryList_brain$n10_short_time$summary_up$enrichr_summary_source$enrichr_5pub$LINCS_L1000_Chem_Pert_Consensus_Sigs)
NeuralSystemCells_down <- filter_results(gene_summaryList_brain$n10_short_time$summary_down$enrichr_summary_source$enrichr_5pub$LINCS_L1000_Chem_Pert_Consensus_Sigs)

BloodCells_up   <- filter_results(gene_summaryList_blood$n10_short_time$summary_up$enrichr_summary_source$enrichr_5pub$LINCS_L1000_Chem_Pert_Consensus_Sigs)
BloodCells_down <- filter_results(gene_summaryList_blood$n10_short_time$summary_down$enrichr_summary_source$enrichr_5pub$LINCS_L1000_Chem_Pert_Consensus_Sigs)

PulmonaryCells_up   <- filter_results(gene_summaryList_lung$n10_short_time$summary_up$enrichr_summary_source$enrichr_5pub$LINCS_L1000_Chem_Pert_Consensus_Sigs)
PulmonaryCells_down <- filter_results(gene_summaryList_lung$n10_short_time$summary_down$enrichr_summary_source$enrichr_5pub$LINCS_L1000_Chem_Pert_Consensus_Sigs)

wb <- createWorkbook()
datasets <- list(
  NeuralSystemCells_up, NeuralSystemCells_down,
  BloodCells_up, BloodCells_down,
  PulmonaryCells_up, PulmonaryCells_down
)

for (i in seq_along(sheet_names)) {
  addWorksheet(wb, sheet_names[i])
  writeData(wb, sheet_names[i], datasets[[i]])
  freezePane(wb, sheet = sheet_names[i], firstRow = TRUE)
  setColWidths(wb, sheet = sheet_names[i], cols = 1:ncol(datasets[[i]]), widths = "auto")
}

output_file <- file.path(output_dir, "LINCSL1000_functionalValidationSignaturesPublication5_fdr0.05overlap3_16.10.2025.xlsx")
saveWorkbook(wb, file = output_file, overwrite = TRUE)
cat("✅ Zapisano:", output_file, "\n")

# ============================================================
# 🧹 Czyszczenie środowiska
# ============================================================

rm(
  NeuralSystemCells_up, NeuralSystemCells_down,
  BloodCells_up, BloodCells_down,
  PulmonaryCells_up, PulmonaryCells_down,
  sheet_names, datasets, wb,
  output_dir, output_file, filter_results
)
gc()  # zwolnienie pamięci

cat("🧹 Środowisko wyczyszczone, wszystko zapisane poprawnie.\n")
