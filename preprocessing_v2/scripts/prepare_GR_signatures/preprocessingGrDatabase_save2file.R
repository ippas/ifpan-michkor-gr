# Ścieżka zapisu
output_path <- "/home/mateusz/projects/ifpan-michkor-gr/results_v2/grDatabase_preprocessing"
file_name <- "GR_gene_database_preprocessing_2026-01-20.xlsx"
full_path <- file.path(output_path, file_name)

# Przygotowanie tabeli
gr_database_df <- papers_data_preprocessing %>%
  filter(
    !simple_tissue %in% c("Other", "placenta")
  ) %>% 
  mutate(
    simple_tissue_regulation = paste0(simple_tissue, "_", regulation)
  ) %>% 
  select(
    source,
    gene_name,
    ensembl_gene_id,
    ensembl_transcript_id,
    refseq_mrna_id,
    hgnc_symbol,
    alias,
    species,
    tissue,
    simple_tissue,
    cell,
    treatment,
    treatment_type,
    dose,
    environment,
    time,
    comparison,
    regulation,
    log2ratio,
    fdr,
    method,
    info
  ) %>% 
  dplyr::rename(
    statistic_method = method,
    detailed_info = info,
    gene_symbol = gene_name
  )

# Utworzenie pliku Excel
wb <- createWorkbook()
addWorksheet(wb, "GR_gene_database")

# Zapis danych
writeData(wb, "GR_gene_database", gr_database_df)

# Zablokowanie pierwszego wiersza (nagłówki)
freezePane(wb, "GR_gene_database", firstRow = TRUE)

# Opcjonalnie: lekka poprawa czytelności
setColWidths(wb, "GR_gene_database", cols = 1:ncol(gr_database_df), widths = "auto")

# Zapis na dysk
saveWorkbook(wb, full_path, overwrite = TRUE)

full_path