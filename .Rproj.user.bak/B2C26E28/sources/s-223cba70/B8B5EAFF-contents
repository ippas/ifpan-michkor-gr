# read data
lite_grSignatures <- gr_genes_signatures_multi_approach_df %>% 
  filter(signature_name %in% c("universal_up", "universal_down", "brain_up", "brain_down")) %>% 
  dplyr::select(c(hgnc_symbol, signature_name)) %>% 
  set_colnames(c("gene_symbol", "signature_name")) %>% 
  mutate(signature_name = case_when(
    signature_name == "universal_up"   ~ "metasignature_up",
    signature_name == "universal_down" ~ "metasignature_down",
    TRUE ~ signature_name
  ))

metasignatureUp_GWASCatalog <- c(
  "CDKN1A", "CRISPLD2", "IL1R2", "FAM107A", "CXCL5",
  "DCN", "ART3", "IL17RA", "PNMT", "ANGPTL4",
  "CD163", "HP", "MAP3K6", "ZBTB16"
)

metasignatureDown_GWASCatalog <- c(
  "IL11", "CHST3", "MME", "IFITM2", "ZBTB16", "KLRC1",
  "CELA2A", "TRPV6", "TOP2A", "CCR7", "FOSL1", "CSF2",
  "BCL3", "TNFRSF11B", "KLRG1", "GDF15", "IL12B"
)

brainUp_GWASCatalog <- c(
  "F13A1", "CYP1B1", "IL1R1", "TNFRSF11A", "S100A11",
  "ERRFI1", "ZBTB16", "HTRA1", "TNS1", "FKBP5",
  "SPSB1", "ANKDD1B", "RAMP2", "IGF2", "PIM3",
  "RHOU", "PER1", "PABPC1L", "FAM107A", "LMOD1",
  "OMG"
)

brainDown_GWASCatalog <- c(
  "DSP", "CSF1", "LPL", "CES2", "RMI2", "TNXB", "EHD4", "IL1R1", "TNFRSF11B",
  "ACAD10", "CTSC", "FAM171A1", "SLCO4A1", "PIEZO2", "FMO2", "LRRC8E", "DENND2C",
  "LAPTM5", "PTPN7", "LAMA1", "RAC2", "PLAU", "CXCL10", "HSD17B14", "CCL18",
  "CORO1A", "ADAP2", "ZNF385C", "WNT7A", "HPRT1", "FCGR3A", "TDRD12"
)

genebass_mentalHealth_skat %>% 
  .$phenocode


genebass_mentalHealth_skat <-read.delim("data/genebass/mentalHealth_Pvalue_SKAT_0.05.tsv.bgz")

# genebass_mentalHealth_burden <-read.delim("data/genebass/mentalHealth_Pvalue_Burden_0.05.tsv.bgz")

grSignaturesLite_association_genebassSKATp05 <- genebass_mentalHealth_skat %>%
  filter(pvalue < 0.05) %>%
  inner_join(lite_grSignatures, by = "gene_symbol") %>%
  dplyr::select(signature_name, everything()) %>%
  dplyr::select(-pvalue_threshold) %>%
  mutate(
    in_GWASCatalog = case_when(
      signature_name == "metasignature_up"   & gene_symbol %in% metasignatureUp_GWASCatalog   ~ 1,
      signature_name == "metasignature_down" & gene_symbol %in% metasignatureDown_GWASCatalog ~ 1,
      signature_name == "brain_up"           & gene_symbol %in% brainUp_GWASCatalog           ~ 1,
      signature_name == "brain_down"         & gene_symbol %in% brainDown_GWASCatalog         ~ 1,
      TRUE ~ 0
    )
  ) %>%
  select(signature_name, gene_id, gene_symbol, in_GWASCatalog, everything()) %>%
  arrange(pvalue) %>%   group_split(signature_name)

names(grSignaturesLite_association_genebassSKATp05) <- c("brain_down", "brain_up", "metasignature_down", "metasignature_up")

# 4. Utwórz workbook i zapisz arkusze
wb <- createWorkbook()

for (sheetname in c("metasignature_up", "metasignature_down", "brain_up", "brain_down")) {
  df <- grSignaturesLite_association_genebassSKATO05[[sheetname]]
  
  addWorksheet(wb, sheetName = sheetname)
  writeData(wb, sheet = sheetname, x = df, withFilter = TRUE)
  freezePane(wb, sheet = sheetname, firstRow = TRUE)
}

# 5. Zapisz workbook
saveWorkbook(wb, file = "data/genebass/grSignaturesLite_association_genebassSKATO05.xlsx", overwrite = TRUE)


# ##############################################################################
# ---- summary grSignatures genebass ----
# ##############################################################################\
summarize_single_df <- function(df) {
  # Upewnij się, że w danych jest kolumna 'phenotype'
  data.frame(
    n_results = nrow(df),
    n_genes = n_distinct(df$gene_symbol),
    n_phenotypes = n_distinct(df$phenotype),
    n_phenotypes_top100 = n_distinct(head(df, 100)$phenotype),
    n_phenotypes_top25 = n_distinct(head(df, 25)$phenotype)
  )
}

bind_rows(
  lapply(grSignaturesLite_association_genebassSKATO05, summarize_single_df),
  .id = "signature_name"
)

